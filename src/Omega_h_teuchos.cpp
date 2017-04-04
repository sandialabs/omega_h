#include "Omega_h_teuchos.hpp"

namespace Omega_h {

void update_var_compare_opts(VarCompareOpts* opts, Teuchos::ParameterList const& pl) {
  if (pl.isType<std::string>("Kind")) {
    auto type_name = pl.get<std::string>("Kind");
    if (type_name == "None") {
      opts->kind = VarCompareOpts::NONE;
    } else if (type_name == "Relative") {
      opts->kind = VarCompareOpts::RELATIVE;
    } else if (type_name == "Absolute") {
      opts->kind = VarCompareOpts::ABSOLUTE;
    } else {
      Omega_h_fail("unknown variable comparison kind \"%s\"\n", type_name.c_str());
    }
  }
  set_if_given(&opts->tolerance, pl, "Tolerance");
  set_if_given(&opts->floor, pl, "Floor");
}

void update_transfer_opts(TransferOpts* opts, Teuchos::ParameterList const& pl) {
  if (pl.isSublist("Fields")) {
    auto& fields_pl = pl.sublist("Fields");
    for (auto it = fields_pl.begin(), end = fields_pl.end(); it != end; ++it) {
      auto field_name = fields_pl.name(it);
      if (fields_pl.isSublist(field_name)) {
        auto& field_pl = fields_pl.sublist(field_name);
        auto type_name = fields_pl.get<std::string>("Type");
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
            update_var_compare_opts(&convergence,
                field_pl.sublist("Diffusion Convergence"));
          }
          opts->integral_diffuse_map[integral_name] = convergence;
        }
        if (opts->type_map[field_name] == OMEGA_H_MOMENTUM_VELOCITY) {
          std::string momentum_name = "momentum";
          if (fields_pl.isType<std::string>("Momentum")) {
            momentum_name = field_pl.get<std::string>("Momentum");
          }
          opts->velocity_momentum_map[field_name] = momentum_name;
          auto convergence = VarCompareOpts::none();
          if (field_pl.isSublist("Diffusion Convergence")) {
            update_var_compare_opts(&convergence,
                field_pl.sublist("Diffusion Convergence"));
          }
          opts->integral_diffuse_map[momentum_name] = convergence;
        }
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
    if (verbosity_name == "Silent") opts->verbosity = SILENT;
    else if (verbosity_name == "Each Adapt") opts->verbosity = EACH_ADAPT;
    else if (verbosity_name == "Each Rebuild") opts->verbosity = EACH_REBUILD;
    else if (verbosity_name == "Extra Stats") opts->verbosity = EXTRA_STATS;
    else Omega_h_fail("unknown verbosity level \"%s\"\n", verbosity_name.c_str());
  }
  set_if_given(&opts->length_histogram_min, pl, "Length Histogram Min");
  set_if_given(&opts->length_histogram_max, pl, "Length Histogram Max");
  set_if_given(&opts->nlength_histogram_bins, pl, "Length Histogram Bin Count");
  set_if_given(&opts->nquality_histogram_bins, pl, "Quality Histogram Bin Count");
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

}
