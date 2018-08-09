#include "Omega_h_teuchos.hpp"

#include "Omega_h_file.hpp"
#include "Omega_h_scatterplot.hpp"

#include <fstream>

namespace Omega_h {

void update_var_compare_opts(VarCompareOpts* opts, Teuchos::ParameterList& pl) {
  if (pl.isType<std::string>("Type")) {
    auto type_name = pl.get<std::string>("Type");
    if (type_name == "None") {
      opts->type = VarCompareOpts::NONE;
    } else if (type_name == "Relative") {
      opts->type = VarCompareOpts::RELATIVE;
    } else if (type_name == "Absolute") {
      opts->type = VarCompareOpts::ABSOLUTE;
    } else {
      Omega_h_fail(
          "unknown variable comparison type \"%s\"\n", type_name.c_str());
    }
  }
  set_if_given(&opts->tolerance, pl, "Tolerance");
  set_if_given(&opts->floor, pl, "Floor");
}

void update_transfer_opts(TransferOpts* opts, Teuchos::ParameterList& pl) {
  if (pl.isSublist("Fields")) {
    auto& fields_pl = pl.sublist("Fields");
    for (auto it = fields_pl.begin(), end = fields_pl.end(); it != end; ++it) {
      auto field_name = fields_pl.name(it);
      if (fields_pl.isSublist(field_name)) {
        auto& field_pl = fields_pl.sublist(field_name);
        auto type_name = field_pl.get<std::string>("Type");
        if (type_name == "Inherit") {
          opts->type_map[field_name] = OMEGA_H_INHERIT;
        } else if (type_name == "Linear Interp") {
          opts->type_map[field_name] = OMEGA_H_LINEAR_INTERP;
        } else if (type_name == "Metric") {
          opts->type_map[field_name] = OMEGA_H_METRIC;
        } else if (type_name == "Density") {
          opts->type_map[field_name] = OMEGA_H_DENSITY;
        } else if (type_name == "Conserve") {
          opts->type_map[field_name] = OMEGA_H_CONSERVE;
        } else if (type_name == "Momentum Velocity") {
          opts->type_map[field_name] = OMEGA_H_MOMENTUM_VELOCITY;
        } else if (type_name == "Pointwise") {
          opts->type_map[field_name] = OMEGA_H_POINTWISE;
        } else {
          Omega_h_fail("unknown transfer type \"%s\"\n", type_name.c_str());
        }
        if (opts->type_map[field_name] == OMEGA_H_CONSERVE) {
          auto integral_name = field_pl.get<std::string>("Integral");
          opts->integral_map[field_name] = integral_name;
          auto convergence = VarCompareOpts::none();
          if (field_pl.isSublist("Diffusion Convergence")) {
            update_var_compare_opts(
                &convergence, field_pl.sublist("Diffusion Convergence"));
          }
          opts->integral_diffuse_map[integral_name] = convergence;
        }
        if (opts->type_map[field_name] == OMEGA_H_MOMENTUM_VELOCITY) {
          std::string momentum_name = "momentum";
          if (field_pl.isType<std::string>("Momentum")) {
            momentum_name = field_pl.get<std::string>("Momentum");
          }
          opts->velocity_momentum_map[field_name] = momentum_name;
          auto density_name = field_pl.get<std::string>("Density");
          opts->velocity_density_map[field_name] = density_name;
          auto convergence = VarCompareOpts::none();
          if (field_pl.isSublist("Diffusion Convergence")) {
            update_var_compare_opts(
                &convergence, field_pl.sublist("Diffusion Convergence"));
          }
          opts->integral_diffuse_map[momentum_name] = convergence;
        }
      } else {
        Omega_h_fail("expected \"%s\" to be a sublist\n", field_name.c_str());
      }
    }
  }
  set_if_given(&opts->should_conserve_size, pl, "Conserve Size");
  set_if_given(&opts->max_size_steps, pl, "Max Size Steps");
  set_if_given(&opts->min_size_step_ratio, pl, "Min Size Step Ratio");
  set_if_given(&opts->max_size_error_ratio, pl, "Max Size Error Ratio");
}

void update_adapt_opts(AdaptOpts* opts, Teuchos::ParameterList& pl) {
  set_if_given(&opts->min_length_desired, pl, "Min Length Desired");
  set_if_given(&opts->max_length_desired, pl, "Max Length Desired");
  set_if_given(&opts->max_length_allowed, pl, "Max Length Allowed");
  set_if_given(&opts->min_quality_allowed, pl, "Min Quality Allowed");
  set_if_given(&opts->min_quality_desired, pl, "Min Quality Desired");
  set_if_given(&opts->nsliver_layers, pl, "Sliver Layer Count");
  if (pl.isType<std::string>("Verbosity")) {
    auto verbosity_name = pl.get<std::string>("Verbosity");
    if (verbosity_name == "Silent")
      opts->verbosity = SILENT;
    else if (verbosity_name == "Each Adapt")
      opts->verbosity = EACH_ADAPT;
    else if (verbosity_name == "Each Rebuild")
      opts->verbosity = EACH_REBUILD;
    else if (verbosity_name == "Extra Stats")
      opts->verbosity = EXTRA_STATS;
    else
      Omega_h_fail("unknown verbosity level \"%s\"\n", verbosity_name.c_str());
  }
  set_if_given(&opts->length_histogram_min, pl, "Length Histogram Min");
  set_if_given(&opts->length_histogram_max, pl, "Length Histogram Max");
  set_if_given(&opts->nlength_histogram_bins, pl, "Length Histogram Bin Count");
  set_if_given(
      &opts->nquality_histogram_bins, pl, "Quality Histogram Bin Count");
#ifdef OMEGA_H_USE_EGADS
  set_if_given(&opts->should_smooth_snap, pl, "Smooth Snap");
  set_if_given(&opts->snap_smooth_tolerance, pl, "Snap Smooth Tolerance");
#endif
  set_if_given(&opts->should_refine, pl, "Refine");
  set_if_given(&opts->should_coarsen, pl, "Coarsen");
  set_if_given(&opts->should_swap, pl, "Swap");
  set_if_given(&opts->should_coarsen_slivers, pl, "Coarsen Slivers");
  if (pl.isSublist("Transfer")) {
    update_transfer_opts(&opts->xfer_opts, pl.sublist("Transfer"));
  }
}

MetricSource get_metric_source(Teuchos::ParameterList& pl) {
  MetricSource source;
  auto type_name = pl.get<std::string>("Type");
  if (type_name == "Constant") {
    source.type = OMEGA_H_CONSTANT;
  } else if (type_name == "Variation") {
    source.type = OMEGA_H_VARIATION;
  } else if (type_name == "Derivative") {
    source.type = OMEGA_H_DERIVATIVE;
  } else if (type_name == "Given") {
    source.type = OMEGA_H_GIVEN;
  } else if (type_name == "Implied") {
    source.type = OMEGA_H_IMPLIED;
  } else if (type_name == "Proximity") {
    source.type = OMEGA_H_PROXIMITY;
  } else if (type_name == "Curvature") {
    source.type = OMEGA_H_CURVATURE;
  } else {
    Omega_h_fail("unknown metric source type \"%s\"\n", type_name.c_str());
  }
  bool should_scale = true;
  set_if_given(&should_scale, pl, "Scale");
  source.scales = should_scale ? OMEGA_H_SCALES : OMEGA_H_ABSOLUTE;
  if (pl.isType<std::string>("Tag Name")) {
    source.tag_name = pl.get<std::string>("Tag Name");
    check_tag_name(source.tag_name);
  }
  source.knob = 1.0;
  set_if_given(&source.knob, pl, "Knob");
  source.isotropy = OMEGA_H_ANISOTROPIC;
  if (pl.isType<std::string>("Isotropy")) {
    auto isotropy_name = pl.get<std::string>("Isotropy");
    if (isotropy_name == "Anisotropic") {
      source.isotropy = OMEGA_H_ANISOTROPIC;
    } else if (isotropy_name == "Preserve Length") {
      source.isotropy = OMEGA_H_ISO_LENGTH;
    } else if (isotropy_name == "Preserve Size") {
      source.isotropy = OMEGA_H_ISO_SIZE;
    } else {
      Omega_h_fail("unknown isotropy type \"%s\"\n", isotropy_name.c_str());
    }
  }
  return source;
}

void update_metric_input(MetricInput* input, Teuchos::ParameterList& pl) {
  if (pl.isSublist("Sources")) {
    auto& sources_pl = pl.sublist("Sources");
    for (auto it = sources_pl.begin(), end = sources_pl.end(); it != end;
         ++it) {
      auto source_name = sources_pl.name(it);
      if (sources_pl.isSublist(source_name)) {
        input->sources.push_back(
            get_metric_source(sources_pl.sublist(source_name)));
      } else {
        Omega_h_fail("expected \"%s\" to be a sublist\n", source_name.c_str());
      }
    }
  }
  set_if_given(&input->should_limit_lengths, pl, "Limit Lengths");
  set_if_given(&input->max_length, pl, "Max Length");
  set_if_given(&input->min_length, pl, "Min Length");
  set_if_given(&input->should_limit_gradation, pl, "Limit Gradation");
  set_if_given(&input->max_gradation_rate, pl, "Max Gradation Rate");
  set_if_given(&input->gradation_convergence_tolerance, pl,
      "Gradation Convergence Tolerance");
  set_if_given(&input->should_limit_element_count, pl, "Limit Element Count");
  set_if_given(&input->max_element_count, pl, "Max Element Count");
  set_if_given(&input->min_element_count, pl, "Min Element Count");
  set_if_given(&input->element_count_over_relaxation, pl,
      "Element Count Over-Relaxation");
  input->nsmoothing_steps = 0;
  set_if_given(&input->nsmoothing_steps, pl, "Smoothing Step Count");
}

Teuchos::RCP<Teuchos::Comm<int>> make_teuchos_comm(CommPtr comm_osh) {
#ifdef OMEGA_H_USE_MPI
  return Teuchos::RCP<Teuchos::Comm<int>>(
      new Teuchos::MpiComm<int>(comm_osh->get_impl()));
#else
  (void)comm_osh;
  return Teuchos::RCP<Teuchos::Comm<int>>(new Teuchos::SerialComm<int>());
#endif
}

void update_parameters_from_file(std::string const& filepath,
    Teuchos::ParameterList* pl, Teuchos::Comm<int> const& comm) {
  if (ends_with(filepath, ".xml")) {
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        filepath, Teuchos::Ptr<Teuchos::ParameterList>(pl), comm);
  } else if (ends_with(filepath, ".yaml")) {
    Teuchos::updateParametersFromYamlFileAndBroadcast(
        filepath, Teuchos::Ptr<Teuchos::ParameterList>(pl), comm, true);
  } else {
    Omega_h_fail(
        "\"%s\" is not a known parameter list format\n", filepath.c_str());
  }
}

void write_parameters(
    std::ostream& stream, Teuchos::ParameterList const& pl, bool is_yaml) {
  std::ios saved_state(nullptr);
  saved_state.copyfmt(stream);
  stream << std::scientific << std::setprecision(17);
  if (is_yaml) {
    Teuchos::writeParameterListToYamlOStream(pl, stream);
  } else {
    Teuchos::writeParameterListToXmlOStream(pl, stream);
  }
  stream.copyfmt(saved_state);
}

void write_parameters(
    std::string const& filepath, Teuchos::ParameterList const& pl) {
  std::ofstream stream(filepath.c_str());
  if (!stream.is_open()) {
    Omega_h_fail("couldn't open \"%s\" to write\n", filepath.c_str());
  }
  write_parameters(stream, pl, ends_with(filepath, ".yaml"));
}

void check_unused(Teuchos::ParameterList const& pl) {
  for (auto it = pl.begin(), end = pl.end(); it != end; ++it) {
    auto& entry = it->second;
    if (entry.isList()) {
      auto& a = entry.getAny(false);
      auto& sl = Teuchos::any_cast<Teuchos::ParameterList>(a);
      check_unused(sl);
    } else if (!entry.isDefault() && !entry.isUsed()) {
      auto name = pl.name();
      name += "->";
      name += it->first;
      Omega_h_fail("Parameter \"%s\" was not used!", name.c_str());
    }
  }
}

void echo_parameters(std::ostream& stream, Teuchos::ParameterList const& pl) {
  Teuchos::ParameterList::PrintOptions opts;
  std::ios saved_state(nullptr);
  saved_state.copyfmt(stream);
  stream << std::scientific << std::setprecision(17);
  opts.showTypes(true);
  opts.showFlags(true);
  opts.showDoc(true);
  pl.print(std::cout, opts);
  stream.copyfmt(saved_state);
}

void update_class_sets(ClassSets* p_sets, Teuchos::ParameterList& pl) {
  ClassSets& sets = *p_sets;
  for (auto it = pl.begin(), end = pl.end(); it != end; ++it) {
    auto set_name = pl.name(it);
    auto pairs = pl.get<Teuchos::TwoDArray<int>>(set_name);
    if (pairs.getNumCols() != 2) {
      Omega_h_fail(
          "Expected \"%s\" to be an array of int pairs\n", set_name.c_str());
    }
    auto npairs = pairs.getNumRows();
    for (decltype(npairs) i = 0; i < npairs; ++i) {
      auto class_dim = Int(pairs(i, 0));
      auto class_id = LO(pairs(i, 1));
      sets[set_name].push_back({class_dim, class_id});
    }
  }
}

static char const* const assoc_param_names[NSET_TYPES] = {
    "Element Sets", "Side Sets", "Node Sets"};

void update_assoc(Assoc* p_assoc, Teuchos::ParameterList& pl) {
  if (pl.isType<std::string>("File")) {
    auto filepath = pl.get<std::string>("File");
    update_from_file(p_assoc, filepath);
  }
  Assoc& assoc = *p_assoc;
  for (Int set_type = 0; set_type < NSET_TYPES; ++set_type) {
    if (pl.isSublist(assoc_param_names[set_type])) {
      auto& set_type_pl = pl.sublist(assoc_param_names[set_type]);
      update_class_sets(&assoc[std::size_t(set_type)], set_type_pl);
    }
  }
}

void update_tag_set(TagSet* p_tags, Int elem_dim, Teuchos::ParameterList& pl) {
  TagSet& tags = *p_tags;
  std::map<std::string, Int> names2dims;
  names2dims["Element"] = elem_dim;
  names2dims["Side"] = elem_dim - 1;
  names2dims["Edge"] = 1;
  names2dims["Vertex"] = 0;
  names2dims["Node"] = 0;
  for (auto pair : names2dims) {
    auto name = pair.first;
    auto dim = pair.second;
    Teuchos::Array<std::string> tag_names;
    Omega_h::set_if_given(&tag_names, pl, name);
    for (auto& tag_name : tag_names) {
      tags[size_t(dim)].insert(tag_name);
    }
  }
}

Int get_ent_dim_by_name(Mesh* mesh, std::string const& name) {
  if (name == "Element")
    return mesh->dim();
  else if (name == "Side")
    return mesh->dim() - 1;
  else if (name == "Node")
    return 0;
  else if (name == "Edge")
    return 1;
  else if (name == "Face")
    return 2;
  else if (name == "Cell")
    return mesh->dim();
  else
    OMEGA_H_NORETURN(-1);
}

template <Int dim>
static void write_scatterplot_dim(Mesh* mesh, Teuchos::ParameterList& pl) {
  auto filepath = pl.get<std::string>("File");
  Int ent_dim = 0;
  if (pl.isType<std::string>("Entity")) {
    auto ent_str = pl.get<std::string>("Entity");
    ent_dim = get_ent_dim_by_name(mesh, ent_str);
  }
  auto tag_name = pl.get<std::string>("Field");
  auto data = mesh->get_array<Real>(ent_dim, tag_name);
  auto origin = zero_vector<dim>();
  if (pl.isType<Teuchos::Array<double>>("Origin")) {
    auto origin_teuchos = pl.get<Teuchos::Array<double>>("Origin");
    for (Int i = 0; i < dim; ++i) origin[i] = origin_teuchos[i];
  }
  std::string separator = "\t";
  if (ends_with(filepath, "csv")) separator = ", ";
  if (pl.isType<Teuchos::Array<double>>("Direction")) {
    Vector<dim> direction;
    auto direction_teuchos = pl.get<Teuchos::Array<double>>("Direction");
    for (Int i = 0; i < dim; ++i) direction[i] = direction_teuchos[i];
    write_linear_scatterplot(
        filepath, mesh, ent_dim, data, direction, origin, separator);
  } else {
    write_radial_scatterplot(filepath, mesh, ent_dim, data, origin, separator);
  }
}

void write_scatterplot(Mesh* mesh, Teuchos::ParameterList& pl) {
  if (mesh->dim() == 3) write_scatterplot_dim<3>(mesh, pl);
  if (mesh->dim() == 2) write_scatterplot_dim<2>(mesh, pl);
  if (mesh->dim() == 1) write_scatterplot_dim<1>(mesh, pl);
}

}  // namespace Omega_h
