#include "Omega_h_compare.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_cmdline.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_owners.hpp"

namespace Omega_h {

VarCompareOpts VarCompareOpts::zero_tolerance() {
  return VarCompareOpts{VarCompareOpts::ABSOLUTE, 0.0, 0.0};
}

VarCompareOpts VarCompareOpts::defaults() {
  return VarCompareOpts{VarCompareOpts::RELATIVE, 1e-6, 0.0};
}

VarCompareOpts VarCompareOpts::none() {
  return VarCompareOpts{VarCompareOpts::NONE, 0.0, 0.0};
}

VarCompareOpts MeshCompareOpts::tag_opts(
    Int dim, std::string const& name) const {
  auto it = tags2opts[dim].find(name);
  if (it == tags2opts[dim].end()) {
    return VarCompareOpts::none();
  }
  return it->second;
}

MeshCompareOpts MeshCompareOpts::init(
    Mesh const* mesh, VarCompareOpts var_opts) {
  MeshCompareOpts opts;
  for (Int dim = 0; dim <= mesh->dim(); ++dim) {
    for (Int i = 0; i < mesh->ntags(dim); ++i) {
      auto tagbase = mesh->get_tag(dim, i);
      opts.tags2opts[dim][tagbase->name()] = var_opts;
    }
  }
  opts.time_step_opts = var_opts;
  return opts;
}

template <typename T>
struct CompareArrays {
  static bool compare(CommPtr comm, Read<T> a, Read<T> b, VarCompareOpts opts,
      Int ncomps, Int dim, bool verbose) {
    if (opts.type == VarCompareOpts::NONE) return true;
    auto this_rank_matches = (a == b);
    if (comm->reduce_and(this_rank_matches)) return true;
    if (!verbose) return false;
    /* if integer arrays are different, we find the first mismatching
       value and print it out as a place for users to start
       tracking down the issue */
    auto h_a = HostRead<T>(a);
    auto h_b = HostRead<T>(b);
    auto global_start = comm->exscan(GO(h_a.size()), OMEGA_H_SUM);
    I32 rank_cand = ArithTraits<I32>::max();
    if (!this_rank_matches) rank_cand = comm->rank();
    // only the first mismatching rank prints
    auto best_rank = comm->allreduce(rank_cand, OMEGA_H_MIN);
    if (comm->rank() == best_rank) {
      for (LO j = 0; j < h_a.size(); ++j) {
        // and it only prints its first mismatching value
        if (h_a[j] != h_b[j]) {
          auto global_comp = global_start + GO(j);
          auto global_ent = global_comp / ncomps;
          auto comp = global_comp % ncomps;
          std::cout << dimensional_singular_name(dim) << ' ' << global_ent
                    << " comp " << comp << " " << I64(h_a[j])
                    << " != " << I64(h_b[j]) << '\n';
          break;
        }
      }
    }
    return false;
  }
};

Real get_real_diff(Real a, Real b, VarCompareOpts opts) {
  if (opts.type == VarCompareOpts::RELATIVE) {
    return rel_diff_with_floor(a, b, opts.floor);
  } else {
    return std::abs(a - b);
  }
}

bool compare_real(Real a, Real b, VarCompareOpts opts) {
  return get_real_diff(a, b, opts) <= opts.tolerance;
}

template <>
struct CompareArrays<Real> {
  static bool compare(CommPtr comm, Read<Real> a, Read<Real> b,
      VarCompareOpts opts, Int ncomps, Int dim, bool verbose) {
    if (opts.type == VarCompareOpts::NONE) return true;
    auto tol = opts.tolerance;
    auto floor = opts.floor;
    if (opts.type == VarCompareOpts::RELATIVE) {
      if (comm->reduce_and(are_close(a, b, tol, floor))) return true;
    } else {
      if (comm->reduce_and(are_close_abs(a, b, tol))) return true;
    }
    if (!verbose) return false;
    /* if floating point arrays are different, we find the value with the
       largest relative difference and print it out for users to determine
       whether this is actually a serious regression
       (and where in the mesh it is most serious)
       or whether tolerances simply need adjusting */
    auto ah = HostRead<Real>(a);
    auto bh = HostRead<Real>(b);
    LO max_i = -1;
    Real max_diff = 0.0;
    for (LO i = 0; i < ah.size(); ++i) {
      auto diff = get_real_diff(ah[i], bh[i], opts);
      if (diff > max_diff) {
        max_i = i;
        max_diff = diff;
      }
    }
    auto global_start = comm->exscan(GO(ah.size()), OMEGA_H_SUM);
    auto global_max_diff = comm->allreduce(max_diff, OMEGA_H_MAX);
    I32 rank_cand = ArithTraits<I32>::max();
    if (max_diff == global_max_diff) rank_cand = comm->rank();
    auto best_rank = comm->allreduce(rank_cand, OMEGA_H_MIN);
    if (comm->rank() == best_rank) {
      auto global_max_i = global_start + max_i;
      auto ent_global = global_max_i / ncomps;
      auto comp = global_max_i % ncomps;
      auto precision_before = std::cout.precision();
      std::ios::fmtflags flags_before(std::cout.flags());
      std::cout << std::scientific << std::setprecision(15);
      std::cout << "max diff at " << dimensional_singular_name(dim) << " "
                << ent_global << ", comp " << comp << ", values " << ah[max_i]
                << " vs " << bh[max_i] << '\n';
      std::cout.flags(flags_before);
      std::cout.precision(precision_before);
    }
    comm->barrier();
    return false;
  }
};

template <typename T>
bool compare_arrays(CommPtr comm, Read<T> a, Read<T> b, VarCompareOpts opts,
    Int ncomps, Int dim, bool verbose) {
  return CompareArrays<T>::compare(comm, a, b, opts, ncomps, dim, verbose);
}

template <typename T>
static bool compare_copy_data(Int dim, Read<T> a_data, Dist a_dist,
    Read<T> b_data, Dist b_dist, Int ncomps, VarCompareOpts opts,
    bool verbose) {
  if (opts.type == VarCompareOpts::NONE) return true;
  auto a_lin_data = reduce_data_to_owners(a_data, a_dist, ncomps);
  auto b_lin_data = reduce_data_to_owners(b_data, b_dist, ncomps);
  OMEGA_H_CHECK(a_lin_data.size() == b_lin_data.size());
  auto comm = a_dist.parent_comm();
  auto ret =
      compare_arrays(comm, a_lin_data, b_lin_data, opts, ncomps, dim, verbose);
  return ret;
}

static Read<GO> get_local_conn(Mesh* mesh, Int dim, Int low_dim) {
  auto h2l = mesh->ask_down(dim, low_dim);
  auto l_globals = mesh->globals(low_dim);
  auto hl2l_globals = unmap(h2l.ab2b, l_globals, 1);
  return hl2l_globals;
}

Omega_h_Comparison compare_meshes(
    Mesh* a, Mesh* b, MeshCompareOpts const& opts, bool verbose, bool full) {
  OMEGA_H_CHECK(a->comm()->size() == b->comm()->size());
  OMEGA_H_CHECK(a->comm()->rank() == b->comm()->rank());
  auto comm = a->comm();
  auto should_print = verbose && (comm->rank() == 0);
  if (a->family() != b->family()) {
    if (should_print) std::cout << "mesh element families differ\n";
    return OMEGA_H_DIFF;
  }
  if (a->dim() != b->dim()) {
    if (should_print) std::cout << "mesh dimensions differ\n";
    return OMEGA_H_DIFF;
  }
  Omega_h_Comparison result = OMEGA_H_SAME;
  for (Int dim = 0; dim <= a->dim(); ++dim) {
    if (a->nglobal_ents(dim) != b->nglobal_ents(dim)) {
      if (should_print) {
        std::cout << "global " << topological_singular_name(a->family(), dim)
                  << " counts differ\n";
      }
      return OMEGA_H_DIFF;
    }
    if (!full && (0 < dim) && (dim < a->dim())) continue;
    auto a_globals = a->globals(dim);
    auto b_globals = b->globals(dim);
    auto a_dist = copies_to_linear_owners(comm, a_globals);
    auto b_dist = copies_to_linear_owners(comm, b_globals);
    if (dim > 0) {
      auto low_dim = ((full) ? (dim - 1) : (VERT));
      auto a_conn = get_local_conn(a, dim, low_dim);
      auto b_conn = get_local_conn(b, dim, low_dim);
      auto deg = element_degree(a->family(), dim, low_dim);
      auto ok = compare_copy_data(dim, a_conn, a_dist, b_conn, b_dist, deg,
          VarCompareOpts::zero_tolerance(), true);
      if (!ok) {
        if (should_print) {
          std::cout << topological_singular_name(a->family(), dim)
                    << " connectivity doesn't match\n";
        }
        result = OMEGA_H_DIFF;
        continue;
      }
    }
    for (Int i = 0; i < a->ntags(dim); ++i) {
      auto tag = a->get_tag(dim, i);
      auto const& name = tag->name();
      if (!b->has_tag(dim, name)) {
        if (should_print) {
          std::cout << topological_singular_name(a->family(), dim) << " tag \""
                    << name << "\" exists in first mesh but not second\n";
        }
        result = OMEGA_H_DIFF;
        continue;
      }
      auto ncomps = tag->ncomps();
      auto tag_opts = opts.tag_opts(dim, name);
      bool ok = false;
      switch (tag->type()) {
        case OMEGA_H_I8:
          ok = compare_copy_data(dim, a->get_array<I8>(dim, name), a_dist,
              b->get_array<I8>(dim, name), b_dist, ncomps, tag_opts, verbose);
          break;
        case OMEGA_H_I32:
          ok = compare_copy_data(dim, a->get_array<I32>(dim, name), a_dist,
              b->get_array<I32>(dim, name), b_dist, ncomps, tag_opts, verbose);
          break;
        case OMEGA_H_I64:
          ok = compare_copy_data(dim, a->get_array<I64>(dim, name), a_dist,
              b->get_array<I64>(dim, name), b_dist, ncomps, tag_opts, verbose);
          break;
        case OMEGA_H_F64:
          ok = compare_copy_data(dim, a->get_array<Real>(dim, name), a_dist,
              b->get_array<Real>(dim, name), b_dist, ncomps, tag_opts, verbose);
          break;
      }
      if (!ok) {
        if (should_print) {
          std::cout << topological_singular_name(a->family(), dim) << " tag \""
                    << name << "\" values are different\n";
        }
        comm->barrier();
        result = OMEGA_H_DIFF;
      }
    }
    for (Int i = 0; i < b->ntags(dim); ++i) {
      auto tag = b->get_tag(dim, i);
      if (!a->has_tag(dim, tag->name())) {
        if (should_print) {
          std::cout << topological_singular_name(a->family(), dim) << " tag \""
                    << tag->name()
                    << "\" exists in second mesh but not in first\n";
        }
        if (result == OMEGA_H_SAME) {
          result = OMEGA_H_MORE;
        }
      }
    }
  }
  return result;
}

void get_diff_program_cmdline(
    std::string const& a_name, std::string const& b_name, CmdLine* p_cmdline) {
  p_cmdline->add_arg<std::string>(a_name);
  p_cmdline->add_arg<std::string>(b_name);
  auto& tolflag = p_cmdline->add_flag(
      "-tolerance", "Overrides the default tolerance of 1.0E-6");
  tolflag.add_arg<double>("value");
  auto& floorflag = p_cmdline->add_flag(
      "-Floor", "Overrides the default floor tolerance of 0.0");
  floorflag.add_arg<double>("value");
  p_cmdline->add_flag("-superset",
      std::string("Allow ") + b_name + " to have more variables than" + a_name);
  auto& fileflag = p_cmdline->add_flag("-f", "Read exodiff command file");
  fileflag.add_arg<std::string>("cmd_file");
}

static VarCompareOpts parse_compare_opts(std::vector<std::string> const& tokens,
    std::size_t start, VarCompareOpts default_opts, std::string const& path,
    Int linenum) {
  auto opts = default_opts;
  if (start == tokens.size()) return opts;
  if (tokens[start] == "relative") {
    opts.type = VarCompareOpts::RELATIVE;
    if (start + 1 == tokens.size()) {
      Omega_h_fail(
          "\"relative\" with no value at %s +%d\n", path.c_str(), linenum);
    }
    opts.tolerance = atof(tokens[start + 1].c_str());
    if (start + 2 < tokens.size() && tokens[start + 2] == "floor") {
      if (start + 3 == tokens.size()) {
        Omega_h_fail(
            "\"floor\" with no value at %s +%d\n", path.c_str(), linenum);
      }
      opts.floor = atof(tokens[start + 3].c_str());
    }
  } else if (tokens[start] == "absolute") {
    opts.type = VarCompareOpts::ABSOLUTE;
    if (start + 1 == tokens.size()) {
      Omega_h_fail(
          "\"absolute\" with no value at %s +%d\n", path.c_str(), linenum);
    }
    opts.tolerance = atof(tokens[start + 1].c_str());
  }
  return opts;
}

static Int parse_dim_token(std::string const& token, Int mesh_dim) {
  if (token == "NODAL") return 0;
  if (token == "ELEMENT") return mesh_dim;
  return -1;
}

static void parse_exodiff_cmd_file(Mesh const* mesh, std::string const& path,
    MeshCompareOpts* p_opts, VarCompareOpts all_defaults) {
  std::ifstream file(path.c_str());
  if (!file.is_open()) {
    Omega_h_fail("Could not open exodiff file \"%s\"\n", path.c_str());
  }
  auto mesh_dim = mesh->dim();
  Int linenum = 1;
  Int variables_dim = -1;
  VarCompareOpts variables_opts;
  for (std::string line; std::getline(file, line);) {
    std::stringstream line_stream(line);
    if (line.empty()) continue;
    std::vector<std::string> tokens;
    for (std::string token; line_stream >> token;) tokens.push_back(token);
    if (tokens.empty()) continue;
    if (tokens[0][0] == '#') continue;
    if (variables_dim != -1 && line[0] == '\t') {
      auto name = tokens[0];
      if (name[0] == '!') {
        name = name.substr(1, std::string::npos);
        if (!mesh->has_tag(variables_dim, name)) {
          Omega_h_fail("variable \"%s\" not in mesh (%s +%d)\n", name.c_str(),
              path.c_str(), linenum);
        }
        auto it = p_opts->tags2opts[variables_dim].find(name);
        if (it == p_opts->tags2opts[variables_dim].end()) {
          Omega_h_fail(
              "directive \"%s\" but %s not included (%s +%d)\n"
              "make sure \"(all)\" is added to the variables line\n",
              tokens[0].c_str(), name.c_str(), path.c_str(), linenum);
        }
        p_opts->tags2opts[variables_dim].erase(it);
      } else {
        if (!mesh->has_tag(variables_dim, name)) {
          Omega_h_fail("variable \"%s\" not in mesh (%s +%d)\n", name.c_str(),
              path.c_str(), linenum);
        }
        auto variable_opts =
            parse_compare_opts(tokens, 1, variables_opts, path, linenum);
        p_opts->tags2opts[variables_dim][name] = variable_opts;
      }
    } else if (variables_dim != -1 && line[0] != '\t' &&
               p_opts->tags2opts[variables_dim].empty()) {
      for (Int i = 0; i < mesh->ntags(variables_dim); ++i) {
        auto tagbase = mesh->get_tag(variables_dim, i);
        p_opts->tags2opts[variables_dim][tagbase->name()] = variables_opts;
      }
      variables_dim = -1;
    } else if (tokens[0] == "COORDINATES") {
      auto opts = parse_compare_opts(tokens, 1, all_defaults, path, linenum);
      p_opts->tags2opts[VERT]["coordinates"] = opts;
    } else if (tokens[0] == "TIME" && tokens.size() > 1 &&
               tokens[1] == "STEPS") {
      auto opts =
          parse_compare_opts(tokens, 2, p_opts->time_step_opts, path, linenum);
      p_opts->time_step_opts = opts;
    } else if (variables_dim == -1 &&
               parse_dim_token(tokens[0], mesh_dim) != -1) {
      variables_dim = parse_dim_token(tokens[0], mesh_dim);
      if (!(tokens.size() > 1 && tokens[1] == "VARIABLES")) {
        Omega_h_fail("bad variables header at %s +%d\n", path.c_str(), linenum);
      }
      std::size_t start = 2;
      bool do_all = false;
      if (tokens.size() > 2 && tokens[2] == "(all)") {
        do_all = true;
        ++start;
      }
      variables_opts =
          parse_compare_opts(tokens, start, all_defaults, path, linenum);
      if (do_all) {
        for (Int i = 0; i < mesh->ntags(variables_dim); ++i) {
          auto tagbase = mesh->get_tag(variables_dim, i);
          p_opts->tags2opts[variables_dim][tagbase->name()] = variables_opts;
        }
      }
    }
    ++linenum;
  }
  if (variables_dim != -1 && p_opts->tags2opts[variables_dim].empty()) {
    for (Int i = 0; i < mesh->ntags(variables_dim); ++i) {
      auto tagbase = mesh->get_tag(variables_dim, i);
      p_opts->tags2opts[variables_dim][tagbase->name()] = variables_opts;
    }
  }
}

void accept_diff_program_cmdline(CmdLine const& cmdline, Mesh const* mesh,
    MeshCompareOpts* p_opts, Omega_h_Comparison* p_max_result) {
  auto all_defaults = VarCompareOpts::defaults();
  if (cmdline.parsed("-tolerance")) {
    all_defaults.tolerance = cmdline.get<double>("-tolerance", "value");
  }
  if (cmdline.parsed("-Floor")) {
    all_defaults.floor = cmdline.get<double>("-Floor", "value");
  }
  if (cmdline.parsed("-superset")) {
    *p_max_result = OMEGA_H_MORE;
  } else {
    *p_max_result = OMEGA_H_SAME;
  }
  if (cmdline.parsed("-f")) {
    parse_exodiff_cmd_file(
        mesh, cmdline.get<std::string>("-f", "cmd_file"), p_opts, all_defaults);
  } else {
    *p_opts = MeshCompareOpts::init(mesh, all_defaults);
  }
}

#define EXPL_INST(T)                                                           \
  template bool compare_arrays(CommPtr comm, Read<T> a, Read<T> b,             \
      VarCompareOpts opts, Int ncomps, Int dim, bool verbose);
EXPL_INST(I8)
EXPL_INST(I32)
EXPL_INST(I64)
EXPL_INST(Real)
#undef EXPL_INST

}  // end namespace Omega_h
