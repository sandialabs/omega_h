#include "Omega_h_teuchos.hpp"

#include "Omega_h_file.hpp"

#include <fstream>

#include <TeuchosCore_config.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#ifdef OMEGA_H_USE_YAML
#include <Teuchos_YamlParameterListHelpers.hpp>
#endif
#include <Teuchos_TwoDArray.hpp>

#ifdef OMEGA_H_USE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif

namespace Omega_h {

void update_var_compare_opts(
    VarCompareOpts* opts, Teuchos::ParameterList const& pl) {
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

void update_transfer_opts(
    TransferOpts* opts, Teuchos::ParameterList const& pl) {
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
}

void update_adapt_opts(AdaptOpts* opts, Teuchos::ParameterList const& pl) {
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
  set_if_given(&opts->max_motion_steps, pl, "Max Motion Steps");
  set_if_given(&opts->motion_step_size, pl, "Motion Step Size");
  set_if_given(&opts->should_refine, pl, "Refine");
  set_if_given(&opts->should_coarsen, pl, "Coarsen");
  set_if_given(&opts->should_swap, pl, "Swap");
  set_if_given(&opts->should_coarsen_slivers, pl, "Coarsen Slivers");
  set_if_given(&opts->should_move_for_quality, pl, "Move For Quality");
  set_if_given(&opts->should_allow_pinching, pl, "Allow Pinching");
  if (pl.isSublist("Transfer")) {
    update_transfer_opts(&opts->xfer_opts, pl.sublist("Transfer"));
  }
}

MetricSource get_metric_source(Teuchos::ParameterList const& pl) {
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

void update_metric_input(MetricInput* input, Teuchos::ParameterList const& pl) {
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
  }
#ifdef OMEGA_H_USE_YAML
  else if (ends_with(filepath, ".yaml")) {
    Teuchos::updateParametersFromYamlFileAndBroadcast(
        filepath, Teuchos::Ptr<Teuchos::ParameterList>(pl), comm);
  }
#endif
  else {
    Omega_h_fail(
        "\"%s\" is not a known parameter list format\n", filepath.c_str());
  }
}

void write_parameters(
    std::string const& filepath, Teuchos::ParameterList const& pl) {
  std::ofstream stream(filepath.c_str());
  if (!stream.is_open()) {
    Omega_h_fail("couldn't open \"%s\" to write\n", filepath.c_str());
  }
  stream << std::scientific << std::setprecision(17);
  if (ends_with(filepath, ".xml")) {
    Teuchos::writeParameterListToXmlOStream(pl, stream);
  }
#ifdef OMEGA_H_USE_YAML
  else if (ends_with(filepath, ".yaml")) {
    Teuchos::writeParameterListToYamlOStream(pl, stream);
  }
#endif
  else {
    Omega_h_fail(
        "\"%s\" is not a known parameter list format\n", filepath.c_str());
  }
}

static char const* const assoc_param_names[NSET_TYPES] = {
    "Element Sets", "Side Sets", "Node Sets"};

void update_assoc(Assoc* p_assoc, Teuchos::ParameterList const& pl) {
  if (pl.isType<std::string>("File")) {
    auto filepath = pl.get<std::string>("File");
    update_from_file(p_assoc, filepath);
  }
  Assoc& assoc = *p_assoc;
  for (Int set_type = 0; set_type < NSET_TYPES; ++set_type) {
    if (pl.isSublist(assoc_param_names[set_type])) {
      auto& set_type_pl = pl.sublist(assoc_param_names[set_type]);
      for (auto it = set_type_pl.begin(), end = set_type_pl.end(); it != end;
           ++it) {
        auto set_name = set_type_pl.name(it);
        auto pairs = set_type_pl.get<Teuchos::TwoDArray<int>>(set_name);
        if (pairs.getNumCols() != 2) {
          Omega_h_fail("Expected \"%s\" to be an array of int pairs\n",
              set_name.c_str());
        }
        auto npairs = pairs.getNumRows();
        for (decltype(npairs) i = 0; i < npairs; ++i) {
          auto class_dim = Int(pairs(i, 0));
          auto class_id = LO(pairs(i, 1));
          assoc[set_type][set_name].push_back({class_dim, class_id});
        }
      }
    }
  }
}

void update_tag_set(
    TagSet* p_tags, Int elem_dim, Teuchos::ParameterList const& pl) {
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

}  // namespace Omega_h
