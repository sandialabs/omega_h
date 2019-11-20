#include "Omega_h_conserve.hpp"

#include <array>
#include <iostream>

#include "Omega_h_adj.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_functors.hpp"
#include "Omega_h_graph.hpp"
#include "Omega_h_host_few.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_r3d.hpp"
#include "Omega_h_transfer.hpp"

namespace Omega_h {

/* we separate sub-cavities even further than by classification: we also
   separate the elements touching the boundary from those that don't. this
   ensures that size error, which can only be created and repaired at the
   boundary, remains at the boundary. it should also be beneficial to do this
   for density errors, because a boundary collapse creates a correspondingly
   large integral error for each size error, and so those integral errors will
   be corrected in the same elements in which the size errors are corrected.
 */

enum BdryStatus {
  NOT_BDRY,  // None of the vertices adjacent to the cavity are on any boundary
  TOUCH_BDRY,  // Adjacent vertices are on the boundary, but the key entity is
               // not
  KEY_BDRY,    // The key entity is on the boundary
};

enum ColorMethod {
  NO_COLOR,          // all elements have the same color
  CLASS_COLOR,       // elements are colored by class ID
  CLASS_BDRY_COLOR,  // by class ID and by whether the element touches the
                     // boundary
};

struct Cavs {
  Graph keys2old_elems;
  Graph keys2new_elems;
  LO size() {
    auto n = keys2old_elems.nnodes();
    OMEGA_H_CHECK(n == keys2new_elems.nnodes());
    return n;
  }
};

struct SeparationResult {
  Cavs separated;
  Cavs remainder;
};

/* note that this does not actually work such that one Graph
   is one color. it works at the per-cavity level, where a
   single sub-cavity is guaranteed to be single-color.
   this helps us speed up operations because we don't actually
   need the outer loop to be over one color, we just need the
   sub-cavities to be single-color as they're processed */
SeparationResult separate_by_color_once(
    Cavs cavs, LOs old_elem_colors, LOs new_elem_colors) {
  auto keys2old = cavs.keys2old_elems;
  auto keys2new = cavs.keys2new_elems;
  auto nkeys = cavs.size();
  auto old_keep_w = Write<I8>(keys2old.nedges());
  auto new_keep_w = Write<I8>(keys2new.nedges(), I8(1));
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto ob = keys2old.a2ab[key];
    auto oe = keys2old.a2ab[key + 1];
    if (ob == oe) {
      return;
    }
    auto old_first = keys2old.ab2b[ob];
    auto color = old_elem_colors[old_first];
    for (auto ko = ob; ko < oe; ++ko) {
      auto elem = keys2old.ab2b[ko];
      old_keep_w[ko] = (old_elem_colors[elem] == color);
    }
    auto nb = keys2new.a2ab[key];
    auto ne = keys2new.a2ab[key + 1];
    for (auto kn = nb; kn < ne; ++kn) {
      auto elem = keys2new.ab2b[kn];
      new_keep_w[kn] = (new_elem_colors[elem] == color);
    }
  };
  parallel_for(nkeys, f, "separate_by_color");
  auto old_keep = Read<I8>(old_keep_w);
  auto new_keep = Read<I8>(new_keep_w);
  auto separated_old = filter_graph(keys2old, old_keep);
  auto separated_new = filter_graph(keys2new, new_keep);
  auto remainder_old = filter_graph(keys2old, invert_marks(old_keep));
  auto remainder_new = filter_graph(keys2new, invert_marks(new_keep));
  OMEGA_H_CHECK(
      remainder_old.nedges() + separated_old.nedges() == keys2old.nedges());
  OMEGA_H_CHECK(
      remainder_new.nedges() + separated_new.nedges() == keys2new.nedges());
  auto separated = Cavs{separated_old, separated_new};
  auto remainder = Cavs{remainder_old, remainder_new};
  return {separated, remainder};
}

using CavsByColor = std::vector<Cavs>;
using CavsByColorMethod = std::array<CavsByColor, 3>;
using CavsByBdryStatus = std::array<CavsByColorMethod, 3>;

static CavsByColor separate_by_color(
    Cavs cavs, LOs old_elem_colors, LOs new_elem_colors) {
  CavsByColor cavs_by_color;
  while (cavs.keys2old_elems.nedges() || cavs.keys2new_elems.nedges()) {
    auto res = separate_by_color_once(cavs, old_elem_colors, new_elem_colors);
    cavs_by_color.push_back(res.separated);
    cavs = res.remainder;
  }
  return cavs_by_color;
}

static LOs get_elem_class_ids(Mesh* mesh) {
  if (mesh->has_tag(mesh->dim(), "class_id")) {
    return mesh->get_array<ClassId>(mesh->dim(), "class_id");
  } else {
    return LOs(mesh->nelems(), 1);
  }
}

static Read<I8> get_elems_are_bdry(Mesh* mesh) {
  auto vert_class_dims = mesh->get_array<I8>(VERT, "class_dim");
  auto verts_are_int = each_eq_to(vert_class_dims, I8(mesh->dim()));
  auto verts_are_bdry = invert_marks(verts_are_int);
  auto elems_are_bdry = mark_up(mesh, VERT, mesh->dim(), verts_are_bdry);
  return elems_are_bdry;
}

static Cavs unmap_cavs(LOs a2b, Cavs c) {
  return {
      unmap_graph(a2b, c.keys2old_elems), unmap_graph(a2b, c.keys2new_elems)};
}

static CavsByBdryStatus separate_cavities(Mesh* old_mesh, Mesh* new_mesh,
    Cavs cavs, Int key_dim, LOs keys2kds, Graph* keys2doms = nullptr) {
  CavsByBdryStatus out;
  auto old_elems_are_bdry_i8 = get_elems_are_bdry(old_mesh);
  auto old_elems_are_bdry = array_cast<LO>(old_elems_are_bdry_i8);
  auto cavs2nbdry_elems =
      graph_reduce(cavs.keys2old_elems, old_elems_are_bdry, 1, OMEGA_H_SUM);
  auto cavs_are_bdry = each_gt(cavs2nbdry_elems, 0);
  auto cavs_arent_bdry = invert_marks(cavs_are_bdry);
  auto int_cavs2cavs = collect_marked(cavs_arent_bdry);
  out[NOT_BDRY][NO_COLOR].push_back(unmap_cavs(int_cavs2cavs, cavs));
  auto kd_class_dims = old_mesh->get_array<I8>(key_dim, "class_dim");
  auto key_class_dims = read(unmap(keys2kds, kd_class_dims, 1));
  auto keys_are_bdry = each_lt(key_class_dims, I8(old_mesh->dim()));
  auto bdry_cavs2cavs = collect_marked(keys_are_bdry);
  out[KEY_BDRY][NO_COLOR].push_back(unmap_cavs(bdry_cavs2cavs, cavs));
  if (keys2doms) *keys2doms = unmap_graph(bdry_cavs2cavs, *keys2doms);
  auto keys_arent_bdry = invert_marks(keys_are_bdry);
  auto cavs_touch_bdry = land_each(cavs_are_bdry, keys_arent_bdry);
  auto touch_cavs2cavs = collect_marked(cavs_touch_bdry);
  out[TOUCH_BDRY][NO_COLOR].push_back(unmap_cavs(touch_cavs2cavs, cavs));
  out[NOT_BDRY][CLASS_COLOR].push_back(out[NOT_BDRY][NO_COLOR][0]);
  out[TOUCH_BDRY][CLASS_COLOR].push_back(out[TOUCH_BDRY][NO_COLOR][0]);
  auto old_elem_class_ids = get_elem_class_ids(old_mesh);
  auto new_elem_class_ids = get_elem_class_ids(new_mesh);
  out[KEY_BDRY][CLASS_COLOR] = separate_by_color(
      out[KEY_BDRY][NO_COLOR][0], old_elem_class_ids, new_elem_class_ids);
  out[NOT_BDRY][CLASS_BDRY_COLOR].push_back(out[NOT_BDRY][NO_COLOR][0]);
  auto old_bdry_colors = old_elems_are_bdry;
  auto new_elems_are_bdry = get_elems_are_bdry(new_mesh);
  auto new_bdry_colors = array_cast<LO>(new_elems_are_bdry);
  out[TOUCH_BDRY][CLASS_BDRY_COLOR] = separate_by_color(
      out[TOUCH_BDRY][NO_COLOR][0], old_bdry_colors, new_bdry_colors);
  out[KEY_BDRY][CLASS_BDRY_COLOR] = separate_by_color(
      out[KEY_BDRY][NO_COLOR][0], old_elem_class_ids, new_elem_class_ids);
  return out;
}

/* this function acts on CLASS_BDRY_COLOR cavities,
   adding in any existing error that was present in a sub-cavity
   and distributes the total to the new sub-cavity */

static void carry_class_bdry_integ_error(Mesh* old_mesh, Mesh* new_mesh,
    Cavs class_bdry_cavs, std::string const& error_name,
    Write<Real> new_elem_errors_w) {
  auto dim = old_mesh->dim();
  auto keys2old = class_bdry_cavs.keys2old_elems;
  auto keys2new = class_bdry_cavs.keys2new_elems;
  auto new_elem_sizes = new_mesh->ask_sizes();
  auto new_cav_elem_sizes = read(unmap(keys2new.ab2b, new_elem_sizes, 1));
  auto old_tag = old_mesh->get_tag<Real>(dim, error_name);
  auto ncomps = old_tag->ncomps();
  auto old_elem_errors = old_tag->array();
  auto old_cav_errors =
      graph_reduce(keys2old, old_elem_errors, ncomps, OMEGA_H_SUM);
  auto new_cav_sizes =
      fan_reduce(keys2new.a2ab, new_cav_elem_sizes, 1, OMEGA_H_SUM);
  auto cav_error_densities =
      divide_each_maybe_zero(old_cav_errors, new_cav_sizes);
  auto cav_elem_error_densities =
      expand(cav_error_densities, keys2new.a2ab, ncomps);
  Reals cav_elem_errors =
      multiply_each(cav_elem_error_densities, new_cav_elem_sizes);
  add_into(cav_elem_errors, keys2new.ab2b, new_elem_errors_w, ncomps);
}

/* given CLASS_COLOR cavities in which we expect the possibility
   of error being introduced, compares the old and new integrals
   and adds the relevant error to the element error array */
static void introduce_class_integ_error(Mesh* old_mesh, Mesh* new_mesh,
    Cavs class_cavs, Reals old_elem_densities, Reals new_elem_densities,
    Write<Real> new_elem_errors_w) {
  auto ncomps =
      divide_no_remainder(old_elem_densities.size(), old_mesh->nelems());
  auto keys2new_elems = class_cavs.keys2new_elems;
  auto keys2old_elems = class_cavs.keys2old_elems;
  auto old_elem_sizes = old_mesh->ask_sizes();
  auto new_elem_sizes = new_mesh->ask_sizes();
  auto old_cav_elem_densities =
      read(unmap(keys2old_elems.ab2b, old_elem_densities, ncomps));
  auto new_cav_elem_densities =
      read(unmap(keys2new_elems.ab2b, new_elem_densities, ncomps));
  auto old_cav_elem_sizes = read(unmap(keys2old_elems.ab2b, old_elem_sizes, 1));
  auto new_cav_elem_sizes = read(unmap(keys2new_elems.ab2b, new_elem_sizes, 1));
  Reals old_cav_elem_integrals =
      multiply_each(old_cav_elem_densities, old_cav_elem_sizes);
  Reals new_cav_elem_integrals =
      multiply_each(new_cav_elem_densities, new_cav_elem_sizes);
  auto old_cav_integrals = fan_reduce(
      keys2old_elems.a2ab, old_cav_elem_integrals, ncomps, OMEGA_H_SUM);
  auto new_cav_integrals = fan_reduce(
      keys2new_elems.a2ab, new_cav_elem_integrals, ncomps, OMEGA_H_SUM);
  auto cav_errors = subtract_each(new_cav_integrals, old_cav_integrals);
  auto new_cav_sizes =
      fan_reduce(keys2new_elems.a2ab, new_cav_elem_sizes, 1, OMEGA_H_SUM);
  auto cav_error_densities = divide_each_maybe_zero(cav_errors, new_cav_sizes);
  auto cav_elem_error_densities =
      expand(cav_error_densities, keys2new_elems.a2ab, ncomps);
  Reals cav_elem_errors =
      multiply_each(cav_elem_error_densities, new_cav_elem_sizes);
  add_into(cav_elem_errors, keys2new_elems.ab2b, new_elem_errors_w, ncomps);
}

struct ConservedBools {
  std::array<bool, 3> this_time;
  std::array<bool, 3> always;
};

static void transfer_integ_error(Mesh* old_mesh, Mesh* new_mesh,
    CavsByBdryStatus const& cavs, Reals old_elem_densities,
    Reals new_elem_densities, std::string const& error_name,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    ConservedBools conserved_bools) {
  begin_code("transfer_integ_error");
  auto dim = old_mesh->dim();
  auto old_tag = old_mesh->get_tag<Real>(dim, error_name);
  auto ncomps = old_tag->ncomps();
  auto old_elem_errors = old_tag->array();
  auto new_elem_errors_w = Write<Real>(new_mesh->nelems() * ncomps, 0.0);
  auto same_errors = read(unmap(same_ents2old_ents, old_elem_errors, ncomps));
  map_into(same_errors, same_ents2new_ents, new_elem_errors_w, ncomps);
  for (std::size_t i = 0; i < 3; ++i) {
    if (!conserved_bools.this_time[i]) {
      for (auto class_cavs : cavs[i][CLASS_COLOR]) {
        introduce_class_integ_error(old_mesh, new_mesh, class_cavs,
            old_elem_densities, new_elem_densities, new_elem_errors_w);
      }
    }
    if (!conserved_bools.always[i]) {
      for (auto bdry_cavs : cavs[i][CLASS_COLOR]) {
        carry_class_bdry_integ_error(
            old_mesh, new_mesh, bdry_cavs, error_name, new_elem_errors_w);
      }
    }
  }
  auto new_elem_errors = Reals(new_elem_errors_w);
  new_mesh->add_tag(dim, error_name, ncomps, new_elem_errors);
  end_code();
}

static void transfer_density_error(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, CavsByBdryStatus const& cavs, TagBase const* tagbase,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    ConservedBools conserved_bools) {
  auto dim = old_mesh->dim();
  auto density_name = tagbase->name();
  if (!opts.integral_map.count(density_name)) {
    Omega_h_fail("integral_map[\"%s\"] undefined\n", density_name.c_str());
  }
  auto integral_name = opts.integral_map.find(density_name)->second;
  auto error_name = integral_name + "_error";
  auto old_tag = as<Real>(tagbase);
  auto old_elem_densities = old_tag->array();
  auto new_elem_densities = new_mesh->get_array<Real>(dim, tagbase->name());
  transfer_integ_error(old_mesh, new_mesh, cavs, old_elem_densities,
      new_elem_densities, error_name, same_ents2old_ents, same_ents2new_ents,
      conserved_bools);
}

static void transfer_momentum_error(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, CavsByBdryStatus const& cavs, TagBase const* tagbase,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    ConservedBools conserved_bools) {
  auto dim = old_mesh->dim();
  auto ncomps = tagbase->ncomps();
  auto velocity_name = tagbase->name();
  auto momentum_name = opts.velocity_momentum_map.find(velocity_name)->second;
  auto density_name = opts.velocity_density_map.find(velocity_name)->second;
  auto error_name = momentum_name + "_error";
  auto old_elem_densities = old_mesh->get_array<Real>(dim, density_name);
  auto new_elem_densities = new_mesh->get_array<Real>(dim, density_name);
  auto old_vert_velocities = old_mesh->get_array<Real>(VERT, velocity_name);
  auto new_vert_velocities = new_mesh->get_array<Real>(VERT, velocity_name);
  auto old_elem_velocities =
      average_field(old_mesh, dim, ncomps, old_vert_velocities);
  auto new_elem_velocities =
      average_field(new_mesh, dim, ncomps, new_vert_velocities);
  Reals old_elem_momenta =
      multiply_each(old_elem_velocities, old_elem_densities);
  Reals new_elem_momenta =
      multiply_each(new_elem_velocities, new_elem_densities);
  transfer_integ_error(old_mesh, new_mesh, cavs, old_elem_momenta,
      new_elem_momenta, error_name, same_ents2old_ents, same_ents2new_ents,
      conserved_bools);
}

struct OpConservation {
  OpConservation() {
    density.always[NOT_BDRY] = true;
    density.always[TOUCH_BDRY] = false;
    density.always[KEY_BDRY] = false;
    momentum.always[NOT_BDRY] = false;
    momentum.always[TOUCH_BDRY] = false;
    momentum.always[KEY_BDRY] = false;
  }
  ConservedBools density;
  ConservedBools momentum;
};

static void transfer_conservation_errors(Mesh* old_mesh,
    TransferOpts const& opts, Mesh* new_mesh, CavsByBdryStatus const& cavs,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    OpConservation op_conservation) {
  auto dim = new_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase)) {
      transfer_density_error(old_mesh, opts, new_mesh, cavs, tagbase,
          same_ents2old_ents, same_ents2new_ents, op_conservation.density);
    }
  }
  for (Int i = 0; i < old_mesh->ntags(VERT); ++i) {
    auto tagbase = old_mesh->get_tag(VERT, i);
    if (is_momentum_velocity(old_mesh, opts, VERT, tagbase)) {
      transfer_momentum_error(old_mesh, opts, new_mesh, cavs, tagbase,
          same_ents2old_ents, same_ents2new_ents, op_conservation.momentum);
    }
  }
}

static Cavs form_initial_cavs(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents) {
  auto dim = new_mesh->dim();
  auto kds2old_elems = old_mesh->ask_graph(key_dim, dim);
  auto keys2old_elems = unmap_graph(keys2kds, kds2old_elems);
  auto keys2new_elems = Graph(keys2prods, prods2new_ents);
  return {keys2old_elems, keys2new_elems};
}

void transfer_conserve_refine(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents) {
  if (!should_conserve_any(old_mesh, opts)) return;
  auto init_cavs = form_initial_cavs(
      old_mesh, new_mesh, EDGE, keys2edges, keys2prods, prods2new_ents);
  auto cavs =
      separate_cavities(old_mesh, new_mesh, init_cavs, EDGE, keys2edges);
  OpConservation op_conservation;
  op_conservation.density.this_time[NOT_BDRY] = true;
  op_conservation.density.this_time[TOUCH_BDRY] = true;
  op_conservation.density.this_time[KEY_BDRY] = true;
  op_conservation.momentum.this_time[NOT_BDRY] = true;
  op_conservation.momentum.this_time[TOUCH_BDRY] = true;
  op_conservation.momentum.this_time[KEY_BDRY] = true;
  transfer_conservation_errors(old_mesh, opts, new_mesh, cavs,
      same_ents2old_ents, same_ents2new_ents, op_conservation);
}

/* intersection-based transfer of density fields.
   note that this is only used in single-material cavities,
   and so it should exactly conserve mass in those cases. */
template <Int dim>
void transfer_by_intersection_dim(Mesh* old_mesh, Mesh* new_mesh,
    TagBase const* tagbase, Cavs cavs, Write<Real> new_data_w) {
  auto keys2old_elems = cavs.keys2old_elems;
  auto keys2new_elems = cavs.keys2new_elems;
  auto ncomps = tagbase->ncomps();
  auto old_tag = as<Real>(tagbase);
  auto old_data = old_tag->array();
  auto old_ev2v = old_mesh->ask_elem_verts();
  auto old_coords = old_mesh->coords();
  auto new_ev2v = new_mesh->ask_elem_verts();
  auto new_coords = new_mesh->coords();
  auto nkeys = cavs.size();
  auto f = OMEGA_H_LAMBDA(LO key) {
    for (auto kne = keys2new_elems.a2ab[key];
         kne < keys2new_elems.a2ab[key + 1]; ++kne) {
      auto new_elem = keys2new_elems.ab2b[kne];
      for (Int comp = 0; comp < ncomps; ++comp) {
        new_data_w[new_elem * ncomps + comp] = 0;
      }
      auto new_verts = gather_verts<dim + 1>(new_ev2v, new_elem);
      auto new_points = gather_vectors<dim + 1, dim>(new_coords, new_verts);
      Real total_intersected_size = 0.0;
      for (auto koe = keys2old_elems.a2ab[key];
           koe < keys2old_elems.a2ab[key + 1]; ++koe) {
        auto old_elem = keys2old_elems.ab2b[koe];
        auto old_verts = gather_verts<dim + 1>(old_ev2v, old_elem);
        auto old_points = gather_vectors<dim + 1, dim>(old_coords, old_verts);
        r3d::Polytope<dim> intersection;
        r3d::intersect_simplices(
            intersection, to_r3d(new_points), to_r3d(old_points));
        auto intersection_size = r3d::measure(intersection);
        for (Int comp = 0; comp < ncomps; ++comp) {
          new_data_w[new_elem * ncomps + comp] +=
              intersection_size * old_data[old_elem * ncomps + comp];
        }
        total_intersected_size += intersection_size;
      }
      for (Int comp = 0; comp < ncomps; ++comp) {
        new_data_w[new_elem * ncomps + comp] /= total_intersected_size;
      }
    }  // end loop over new elements
  };
  parallel_for(nkeys, f, "transfer_by_intersection");
}

static void transfer_by_intersection(Mesh* old_mesh, Mesh* new_mesh,
    TagBase const* tagbase, Cavs cavs, Write<Real> new_elem_densities_w) {
  auto dim = old_mesh->dim();
  if (dim == 3) {
    transfer_by_intersection_dim<3>(
        old_mesh, new_mesh, tagbase, cavs, new_elem_densities_w);
  } else if (dim == 2) {
    transfer_by_intersection_dim<2>(
        old_mesh, new_mesh, tagbase, cavs, new_elem_densities_w);
  } else if (dim == 1) {
    transfer_by_intersection_dim<1>(
        old_mesh, new_mesh, tagbase, cavs, new_elem_densities_w);
  } else {
    Omega_h_fail("unsupported dim %d\n", dim);
  }
}

void transfer_densities_and_conserve_swap(Mesh* old_mesh,
    TransferOpts const& opts, Mesh* new_mesh, LOs keys2edges, LOs keys2prods,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents) {
  if (!has_densities_or_conserved(old_mesh, opts)) return;
  auto init_cavs = form_initial_cavs(
      old_mesh, new_mesh, EDGE, keys2edges, keys2prods, prods2new_ents);
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase) ||
        is_density(old_mesh, opts, dim, tagbase)) {
      auto ncomps = tagbase->ncomps();
      auto new_elem_densities_w = Write<Real>(new_mesh->nelems() * ncomps);
      transfer_by_intersection(
          old_mesh, new_mesh, tagbase, init_cavs, new_elem_densities_w);
      transfer_common2(old_mesh, new_mesh, dim, same_ents2old_ents,
          same_ents2new_ents, tagbase, new_elem_densities_w);
    }
  }
  if (!should_conserve_any(old_mesh, opts)) return;
  OpConservation op_conservation;
  op_conservation.density.this_time[NOT_BDRY] = true;
  op_conservation.density.this_time[TOUCH_BDRY] = true;
  op_conservation.density.this_time[KEY_BDRY] = false;
  op_conservation.momentum.this_time[NOT_BDRY] = false;
  op_conservation.momentum.this_time[TOUCH_BDRY] = false;
  op_conservation.momentum.this_time[KEY_BDRY] = false;
  auto cavs =
      separate_cavities(old_mesh, new_mesh, init_cavs, EDGE, keys2edges);
  transfer_conservation_errors(old_mesh, opts, new_mesh, cavs,
      same_ents2old_ents, same_ents2new_ents, op_conservation);
}

void transfer_densities_and_conserve_coarsen(Mesh* old_mesh,
    TransferOpts const& opts, Mesh* new_mesh, LOs keys2verts, Adj keys2doms,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents) {
  if (!has_densities_or_conserved(old_mesh, opts)) return;
  auto keys2prods = keys2doms.a2ab;
  auto init_cavs = form_initial_cavs(
      old_mesh, new_mesh, VERT, keys2verts, keys2prods, prods2new_ents);
  auto bdry_keys2doms = keys2doms;
  auto cavs = separate_cavities(
      old_mesh, new_mesh, init_cavs, VERT, keys2verts, &bdry_keys2doms);
  auto dim = old_mesh->dim();
  for (Int i = 0; i < old_mesh->ntags(dim); ++i) {
    auto tagbase = old_mesh->get_tag(dim, i);
    if (should_conserve(old_mesh, opts, dim, tagbase) ||
        is_density(old_mesh, opts, dim, tagbase)) {
      auto ncomps = tagbase->ncomps();
      auto new_elem_densities_w = Write<Real>(new_mesh->nelems() * ncomps);
      transfer_by_intersection(old_mesh, new_mesh, tagbase,
          cavs[NOT_BDRY][NO_COLOR][0], new_elem_densities_w);
      transfer_by_intersection(old_mesh, new_mesh, tagbase,
          cavs[TOUCH_BDRY][NO_COLOR][0], new_elem_densities_w);
      for (auto color_cavs : cavs[KEY_BDRY][CLASS_COLOR]) {
        transfer_by_intersection(
            old_mesh, new_mesh, tagbase, color_cavs, new_elem_densities_w);
      }
      transfer_common2(old_mesh, new_mesh, dim, same_ents2old_ents,
          same_ents2new_ents, tagbase, new_elem_densities_w);
    }
  }
  if (!should_conserve_any(old_mesh, opts)) return;
  OpConservation op_conservation;
  op_conservation.density.this_time[NOT_BDRY] = true;
  op_conservation.density.this_time[TOUCH_BDRY] = true;
  op_conservation.density.this_time[KEY_BDRY] = false;
  op_conservation.momentum.this_time[NOT_BDRY] = false;
  op_conservation.momentum.this_time[TOUCH_BDRY] = false;
  op_conservation.momentum.this_time[KEY_BDRY] = false;
  transfer_conservation_errors(old_mesh, opts, new_mesh, cavs,
      same_ents2old_ents, same_ents2new_ents, op_conservation);
}

void transfer_conserve_motion(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2verts, Graph keys2elems, LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  if (!should_conserve_any(old_mesh, opts)) return;
  auto keys2prods = keys2elems.a2ab;
  auto prods2new_ents = keys2elems.ab2b;
  auto init_cavs = form_initial_cavs(
      old_mesh, new_mesh, VERT, keys2verts, keys2prods, prods2new_ents);
  OMEGA_H_CHECK(init_cavs.keys2old_elems.a2ab == init_cavs.keys2new_elems.a2ab);
  OMEGA_H_CHECK(init_cavs.keys2old_elems.ab2b == init_cavs.keys2new_elems.ab2b);
  auto bdry_keys2doms = keys2elems;
  auto cavs = separate_cavities(
      old_mesh, new_mesh, init_cavs, VERT, keys2verts, &bdry_keys2doms);
  OpConservation op_conservation;
  op_conservation.density.this_time[NOT_BDRY] = true;
  op_conservation.density.this_time[TOUCH_BDRY] = true;
  op_conservation.density.this_time[KEY_BDRY] = false;
  op_conservation.momentum.this_time[NOT_BDRY] = false;
  op_conservation.momentum.this_time[TOUCH_BDRY] = false;
  op_conservation.momentum.this_time[KEY_BDRY] = false;
  transfer_conservation_errors(old_mesh, opts, new_mesh, cavs,
      same_ents2old_ents, same_ents2new_ents, op_conservation);
}

void fix_momentum_velocity_verts(
    Mesh* mesh, std::vector<ClassPair> const& class_pairs, Int comp) {
  for (Int ent_dim = VERT; ent_dim <= mesh->dim(); ++ent_dim) {
    auto ent_marks = mark_class_closures(mesh, ent_dim, class_pairs);
    auto comp_marks = multiply_each_by(ent_marks, I8(1 << comp));
    if (mesh->has_tag(ent_dim, "momentum_velocity_fixed")) {
      auto old_marks = mesh->get_array<I8>(ent_dim, "momentum_velocity_fixed");
      auto new_marks = bit_or_each(old_marks, comp_marks);
      mesh->set_tag(ent_dim, "momentum_velocity_fixed", new_marks);
    } else {
      mesh->add_tag(ent_dim, "momentum_velocity_fixed", 1, comp_marks);
    }
  }
}

static Read<I8> get_comps_are_fixed(Mesh* mesh) {
  if (mesh->has_tag(VERT, "momentum_velocity_fixed")) {
    return mesh->get_array<I8>(VERT, "momentum_velocity_fixed");
  } else {
    return Read<I8>(mesh->nverts(), I8(0));
  }
}

void setup_conservation_tags(Mesh* mesh, AdaptOpts const& opts) {
  auto& xfer_opts = opts.xfer_opts;
  auto dim = mesh->dim();
  for (Int tagi = 0; tagi < mesh->ntags(dim); ++tagi) {
    auto tagbase = mesh->get_tag(dim, tagi);
    if (should_conserve(mesh, xfer_opts, dim, tagbase)) {
      auto density_name = tagbase->name();
      auto it = xfer_opts.integral_map.find(density_name);
      if (it == xfer_opts.integral_map.end()) {
        Omega_h_fail("conserved density \"%s\" has no integral_map entry\n",
            density_name.c_str());
      }
      auto integral_name = it->second;
      auto error_name = integral_name + "_error";
      auto ncomps = tagbase->ncomps();
      mesh->add_tag(
          dim, error_name, ncomps, Reals(mesh->nelems() * ncomps, 0.0));
    }
  }
  for (Int tagi = 0; tagi < mesh->ntags(VERT); ++tagi) {
    auto tagbase = mesh->get_tag(VERT, tagi);
    if (is_momentum_velocity(mesh, xfer_opts, VERT, tagbase)) {
      auto velocity_name = tagbase->name();
      auto momentum_name =
          xfer_opts.velocity_momentum_map.find(velocity_name)->second;
      auto error_name = momentum_name + "_error";
      auto ncomps = tagbase->ncomps();
      mesh->add_tag(
          dim, error_name, ncomps, Reals(mesh->nelems() * ncomps, 0.0));
    }
  }
}

/* conservative diffusion:
   1) diffusion is not allowed across classification boundaries.
      this preserves conservation on a per-object basis
   2) we are careful with the mathematics of the diffusion algorithm,
      such that it should preserve the total error, and hence
      conservation
 */
static Graph get_elem_diffusion_graph(Mesh* mesh) {
  auto dim = mesh->dim();
  OMEGA_H_CHECK(mesh->owners_have_all_upward(dim - 1));
  auto side_class_dims = mesh->get_array<I8>(dim - 1, "class_dim");
  auto sides_are_bdry = each_lt(side_class_dims, I8(dim));
  if (mesh->comm()->size() > 1) {
    /* in a parallel mesh we need to mark the mesh sides on partition
       boundaries as stopping points as well, otherwise graph construction
       will try to reach across them.
       ghosted elements will make redudant, erroneous computations
       that will be overwritten by their owned copy */
    sides_are_bdry = lor_each(sides_are_bdry, mark_exposed_sides(mesh));
  }
  auto elems2sides = mesh->ask_down(dim, dim - 1);
  auto sides2elems = mesh->ask_up(dim - 1, dim);
  return elements_across_sides(dim, elems2sides, sides2elems, sides_are_bdry);
}

Reals diffuse_densities_once(
    Mesh* mesh, Graph g, Reals densities, Reals cell_sizes) {
  auto out = deep_copy(densities);
  auto max_deg = mesh->dim() + 1;
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto s = cell_sizes[e];
    for (auto ee = g.a2ab[e]; ee < g.a2ab[e + 1]; ++ee) {
      auto oe = g.ab2b[ee];
      auto os = cell_sizes[oe];
      auto mins = min2(s, os);  // minimum of this and other cell size
      /* get difference in densities, multiply by (mins / max_deg)
         to get a mass value that is below the stable limit,
         then divide that mass value by (s) to get the density
         delta to add to this density */
      auto factor = mins / (s * max_deg);
      auto x = densities[e];
      auto ox = densities[oe];
      out[e] += (ox - x) * factor;
    }
  };
  parallel_for(mesh->nelems(), f, "diffuse_densities");
  return mesh->sync_array(mesh->dim(), Reals(out), 1);
}

struct AllBounded : public AndFunctor {
  Reals a;
  Real b;
  AllBounded(Reals a_, Real b_) : a(a_), b(b_) {}
  OMEGA_H_DEVICE void operator()(LO i, value_type& update) const {
    update = update && (std::abs(a[i]) <= b);
  }
};

static bool all_bounded(CommPtr comm, Reals a, Real b) {
  return bool(get_min(comm, each_leq_to(fabs_each(a), b)));
}

static Reals diffuse_densities(Mesh* mesh, Graph g, Reals densities,
    Reals cell_sizes, VarCompareOpts opts, std::string const& name,
    bool verbose) {
  auto comm = mesh->comm();
  Int niters = 0;
  for (niters = 0; !all_bounded(comm, densities, opts.tolerance); ++niters) {
    densities = diffuse_densities_once(mesh, g, densities, cell_sizes);
  }
  if (verbose && !comm->rank()) {
    std::cout << "diffused " << name << " in " << niters << " iterations\n";
  }
  return densities;
}

static Reals diffuse_integrals_weighted(Mesh* mesh, Graph g,
    Reals error_integrals, Reals quantity_integrals, VarCompareOpts opts,
    std::string const& name, bool verbose) {
  if (opts.type == VarCompareOpts::NONE) return error_integrals;
  auto ncomps = divide_no_remainder(error_integrals.size(), g.nnodes());
  if (ncomps > 1) {
    Write<Real> out(error_integrals.size());
    for (Int c = 0; c < ncomps; ++c) {
      auto comp_integrals = get_component(error_integrals, ncomps, c);
      auto comp_quantity_integrals =
          get_component(quantity_integrals, ncomps, c);
      auto comp_name = name + "_" + std::to_string(c);
      comp_integrals = diffuse_integrals_weighted(mesh, g, comp_integrals,
          comp_quantity_integrals, opts, comp_name, verbose);
      set_component(out, comp_integrals, ncomps, c);
    }
    return out;
  }
  auto weighted_sizes = fabs_each(quantity_integrals);
  weighted_sizes = each_max_with(weighted_sizes, opts.floor);
  auto weighted_densities =
      divide_each_maybe_zero(error_integrals, weighted_sizes);
  weighted_densities = diffuse_densities(
      mesh, g, weighted_densities, weighted_sizes, opts, name, verbose);
  error_integrals = multiply_each(weighted_densities, weighted_sizes);
  return error_integrals;
}

static void correct_density_error(Mesh* mesh, TransferOpts const& xfer_opts,
    Graph diffusion_graph, TagBase const* tagbase, bool verbose) {
  auto dim = mesh->dim();
  auto density_name = tagbase->name();
  auto integral_name = xfer_opts.integral_map.find(density_name)->second;
  auto error_name = integral_name + "_error";
  auto ncomps = tagbase->ncomps();
  auto old_densities = mesh->get_array<Real>(dim, density_name);
  mesh->add_tag(dim, std::string("old_") + density_name, ncomps, old_densities);
  auto sizes = mesh->ask_sizes();
  Reals old_integrals = multiply_each(old_densities, sizes);
  auto errors = mesh->get_array<Real>(dim, error_name);
  auto diffuse_tol = xfer_opts.integral_diffuse_map.find(integral_name)->second;
  errors = diffuse_integrals_weighted(mesh, diffusion_graph, errors,
      old_integrals, diffuse_tol, error_name, verbose);
  mesh->set_tag(dim, error_name, errors);
  auto new_integrals = subtract_each(old_integrals, errors);
  auto new_densities = read(divide_each(new_integrals, sizes));
  mesh->set_tag(dim, density_name, new_densities);
  mesh->remove_tag(dim, error_name);
}

void correct_momentum_error(Mesh* mesh, TransferOpts const& xfer_opts,
    Graph diffusion_graph, TagBase const* tagbase, bool verbose) {
  auto dim = mesh->dim();
  auto ncomps = tagbase->ncomps();
  auto velocity_name = tagbase->name();
  auto momentum_name =
      xfer_opts.velocity_momentum_map.find(velocity_name)->second;
  auto density_name =
      xfer_opts.velocity_density_map.find(velocity_name)->second;
  auto error_name = momentum_name + "_error";
  auto elem_densities = mesh->get_array<Real>(dim, density_name);
  auto elem_sizes = mesh->ask_sizes();
  Reals elem_masses = multiply_each(elem_densities, elem_sizes);
  auto vert_velocities = mesh->get_array<Real>(VERT, velocity_name);
  auto old_elem_densities =
      mesh->get_array<Real>(dim, std::string("old_") + density_name);
  Reals old_elem_masses = multiply_each(old_elem_densities, elem_sizes);
  auto elem_velocities = average_field(mesh, dim, ncomps, vert_velocities);
  Reals old_elem_momenta = multiply_each(elem_velocities, old_elem_masses);
  Reals new_elem_momenta = multiply_each(elem_velocities, elem_masses);
  auto elem_errors_from_density =
      subtract_each(new_elem_momenta, old_elem_momenta);
  auto verts2elems = mesh->ask_up(VERT, dim);
  auto vert_masses = graph_reduce(verts2elems, elem_masses, 1, OMEGA_H_SUM);
  vert_masses = divide_each_by(vert_masses, Real(dim + 1));
  auto elems2verts = mesh->ask_down(dim, VERT);
  auto all_flags = get_comps_are_fixed(mesh);
  auto elem_errors = mesh->get_array<Real>(dim, error_name);
  elem_errors = add_each(elem_errors, elem_errors_from_density);
  auto diffuse_tol = xfer_opts.integral_diffuse_map.find(momentum_name)->second;
  elem_errors = diffuse_integrals_weighted(mesh, diffusion_graph, elem_errors,
      new_elem_momenta, diffuse_tol, error_name, verbose);
  mesh->set_tag(dim, error_name, elem_errors);
  auto out = deep_copy(vert_velocities);
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto v_flags = all_flags[v];
    auto v_mass = vert_masses[v];
    for (auto ve = verts2elems.a2ab[v]; ve < verts2elems.a2ab[v + 1]; ++ve) {
      auto e = verts2elems.ab2b[ve];
      auto nfree_verts = zero_vector<3>();
      for (auto ev = e * (dim + 1); ev < (e + 1) * (dim + 1); ++ev) {
        auto v2 = elems2verts.ab2b[ev];
        auto v2_flags = all_flags[v2];
        for (Int comp = 0; comp < ncomps; ++comp) {
          if (!(v2_flags & (1 << comp))) nfree_verts[comp] += 1.0;
        }
      }
      for (Int comp = 0; comp < ncomps; ++comp) {
        if (!(v_flags & (1 << comp))) {
          out[v * ncomps + comp] -=
              elem_errors[e * ncomps + comp] / (nfree_verts[comp] * v_mass);
        }
      }
    }
  };
  parallel_for(mesh->nverts(), f, "correct_momentum_error");
  auto new_velocities = Reals(out);
  new_velocities = mesh->sync_array(VERT, new_velocities, ncomps);
  mesh->set_tag(VERT, velocity_name, new_velocities);
  mesh->remove_tag(dim, error_name);
}

void correct_integral_errors(Mesh* mesh, AdaptOpts const& opts) {
  auto& xfer_opts = opts.xfer_opts;
  if (!should_conserve_any(mesh, xfer_opts)) return;
  begin_code("correct_integral_errors");
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto verbose = opts.verbosity > SILENT;
  auto diffusion_graph = get_elem_diffusion_graph(mesh);
  auto elem_sizes = mesh->ask_sizes();
  auto dim = mesh->dim();
  for (Int tagi = 0; tagi < mesh->ntags(dim); ++tagi) {
    auto tagbase = mesh->get_tag(dim, tagi);
    if (should_conserve(mesh, xfer_opts, dim, tagbase)) {
      correct_density_error(mesh, xfer_opts, diffusion_graph, tagbase, verbose);
    }
  }
  for (Int tagi = 0; tagi < mesh->ntags(VERT); ++tagi) {
    auto tagbase = mesh->get_tag(VERT, tagi);
    if (is_momentum_velocity(mesh, xfer_opts, VERT, tagbase)) {
      correct_momentum_error(
          mesh, xfer_opts, diffusion_graph, tagbase, verbose);
    }
  }
  for (Int tagi = 0; tagi < mesh->ntags(dim); ++tagi) {
    auto tagbase = mesh->get_tag(dim, tagi);
    if (should_conserve(mesh, xfer_opts, dim, tagbase)) {
      auto density_name = tagbase->name();
      mesh->remove_tag(dim, std::string("old_") + density_name);
    }
  }
  end_code();
}

}  // end namespace Omega_h
