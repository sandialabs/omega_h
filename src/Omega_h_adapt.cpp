#include <iomanip>
#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_coarsen.hpp"
#include "Omega_h_confined.hpp"
#include "Omega_h_conserve.hpp"
#include "Omega_h_histogram.hpp"
#include "Omega_h_laplace.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_profile.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_refine.hpp"
#include "Omega_h_swap.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_transfer.hpp"

#ifdef OMEGA_H_USE_EGADS
#include "Omega_h_egads.hpp"
#endif

namespace Omega_h {

void UserTransfer::out_of_line_virtual_method() {}

TransferOpts::TransferOpts() {}

void TransferOpts::validate(Mesh* mesh) const {
  for (auto& pair : type_map) {
    auto& name = pair.first;
    bool tag_exists = false;
    for (Int d = 0; d <= mesh->dim(); ++d) {
      if (mesh->has_tag(d, name)) tag_exists = true;
    }
    if (!tag_exists) {
      Omega_h_fail("Field \"%s\" needs to be transferred but is not attached\n",
          name.c_str());
    }
    if (pair.second == OMEGA_H_MOMENTUM_VELOCITY) {
      auto velocity_name = pair.first;
      OMEGA_H_CHECK(velocity_momentum_map.count(velocity_name));
      OMEGA_H_CHECK(velocity_density_map.count(velocity_name));
      auto density_name = velocity_density_map.find(velocity_name)->second;
      OMEGA_H_CHECK(mesh->has_tag(mesh->dim(), density_name));
      auto density = mesh->get_tagbase(mesh->dim(), density_name);
      OMEGA_H_CHECK(density->type() == OMEGA_H_REAL);
      OMEGA_H_CHECK(density->ncomps() == 1);
    }
  }
}

AdaptOpts::AdaptOpts(Int dim) {
  min_length_desired = 1.0 / std::sqrt(2.0);
  max_length_desired = std::sqrt(2.0);
  max_length_allowed = max_length_desired * 2.0;
  if (dim == 3) {
    min_quality_allowed = 0.20;
    min_quality_desired = 0.30;
  } else if (dim == 2) {
    min_quality_allowed = 0.30;
    min_quality_desired = 0.40;
  } else if (dim == 1) {
    min_quality_allowed = 0.0;
    min_quality_desired = 0.0;
  } else {
    Omega_h_fail("unexpected dim %d\n", dim);
  }
  nsliver_layers = 4;
  verbosity = EACH_REBUILD;
  length_histogram_min = 0.0;
  length_histogram_max = 3.0;
  nlength_histogram_bins = 10;
  nquality_histogram_bins = 10;
#ifdef OMEGA_H_USE_EGADS
  egads_model = nullptr;
  should_smooth_snap = true;
  snap_smooth_tolerance = 1e-2;
  allow_snap_failure = false;
#endif
  should_refine = true;
  should_coarsen = true;
  should_swap = true;
  should_coarsen_slivers = true;
  should_prevent_coarsen_flip = false;
}

static Reals get_fixable_qualities(Mesh* mesh, AdaptOpts const&) {
  /* This used to be an attempt to continue adapting when certain
     elements were constrained to by geometry to have small dihedral angles.
     We'll leave it here as a placeholder for reimplementing such a system
     in the future, but for now it just returns all qualities */
  return mesh->ask_qualities();
}

Real min_fixable_quality(Mesh* mesh, AdaptOpts const& opts) {
  return get_min(mesh->comm(), get_fixable_qualities(mesh, opts));
}

AdaptOpts::AdaptOpts(Mesh* mesh) : AdaptOpts(mesh->dim()) {}

static void adapt_summary(Mesh* mesh, AdaptOpts const& opts,
    MinMax<Real> qualstats, MinMax<Real> lenstats) {
  print_goal_stats(mesh, "quality", mesh->dim(),
      get_fixable_qualities(mesh, opts),
      {opts.min_quality_allowed, opts.min_quality_desired}, qualstats);
  print_goal_stats(mesh, "length", EDGE, mesh->ask_lengths(),
      {opts.min_length_desired, opts.max_length_desired}, lenstats);
}

bool print_adapt_status(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_TIME_FUNCTION;
  auto qualstats = get_minmax(mesh->comm(), get_fixable_qualities(mesh, opts));
  auto lenstats = get_minmax(mesh->comm(), mesh->ask_lengths());
  if (opts.verbosity > SILENT) {
    adapt_summary(mesh, opts, qualstats, lenstats);
  }
  return (qualstats.min >= opts.min_quality_desired &&
          lenstats.min >= opts.min_length_desired &&
          lenstats.max <= opts.max_length_desired);
}

void print_adapt_histograms(Mesh* mesh, AdaptOpts const& opts) {
  auto qh = get_histogram(mesh, mesh->dim(), opts.nquality_histogram_bins, 0.0,
      1.0, mesh->ask_qualities());
  auto lh = get_histogram(mesh, EDGE, opts.nlength_histogram_bins,
      opts.length_histogram_min, opts.length_histogram_max,
      mesh->ask_lengths());
  auto owned_qualities =
      mesh->owned_array(mesh->dim(), mesh->ask_qualities(), 1);
  auto qual_sum = get_sum(mesh->comm(), owned_qualities);
  auto global_nelems = mesh->nglobal_ents(mesh->dim());
  auto avg_qual = qual_sum / global_nelems;
  if (can_print(mesh)) {
    print_histogram(qh, "quality");
    print_histogram(lh, "length");
    std::cout << "average quality: " << avg_qual << '\n';
  }
}

static void validate(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_CHECK(0.0 <= opts.min_quality_allowed);
  OMEGA_H_CHECK(opts.min_quality_allowed <= opts.min_quality_desired);
  OMEGA_H_CHECK(opts.min_quality_desired <= 1.0);
  OMEGA_H_CHECK(opts.nsliver_layers >= 0);
  OMEGA_H_CHECK(opts.nsliver_layers < 100);
  auto mq = min_fixable_quality(mesh, opts);
  if (mq < opts.min_quality_allowed && !mesh->comm()->rank()) {
    std::cout << "WARNING: worst input element has quality " << mq
              << " but minimum allowed is " << opts.min_quality_allowed << "\n";
  }
}

static bool pre_adapt(Mesh* mesh, AdaptOpts const& opts) {
  validate(mesh, opts);
  opts.xfer_opts.validate(mesh);
  if (opts.verbosity >= EACH_ADAPT && !mesh->comm()->rank()) {
    std::cout << "before adapting:\n";
  }
  if (print_adapt_status(mesh, opts)) return false;
  if (opts.verbosity >= EXTRA_STATS) {
    print_adapt_histograms(mesh, opts);
  }
  if ((opts.verbosity >= EACH_REBUILD) && !mesh->comm()->rank()) {
    std::cout << "addressing edge lengths\n";
  }
  return true;
}

static void post_rebuild(Mesh* mesh, AdaptOpts const& opts) {
  if (opts.verbosity >= EACH_REBUILD) print_adapt_status(mesh, opts);
}

static void satisfy_lengths(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_TIME_FUNCTION;
  bool did_anything;
  do {
    did_anything = false;
    if (opts.should_refine && refine_by_size(mesh, opts)) {
      post_rebuild(mesh, opts);
      did_anything = true;
    }
    if (opts.should_coarsen && coarsen_by_size(mesh, opts)) {
      post_rebuild(mesh, opts);
      did_anything = true;
    }
  } while (did_anything);
}

static bool satisfy_quality(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_TIME_FUNCTION;
  if (min_fixable_quality(mesh, opts) >= opts.min_quality_desired) return true;
  if ((opts.verbosity >= EACH_REBUILD) && can_print(mesh)) {
    std::cout << "addressing element qualities\n";
  }
  do {
    if (opts.should_swap && swap_edges(mesh, opts)) {
      post_rebuild(mesh, opts);
      continue;
    }
    if (opts.should_coarsen_slivers && coarsen_slivers(mesh, opts)) {
      post_rebuild(mesh, opts);
      continue;
    }
    if ((opts.verbosity > SILENT) && can_print(mesh)) {
      std::cout << "could not satisfy quality\n";
    }
    return false;
  } while (min_fixable_quality(mesh, opts) < opts.min_quality_desired);
  return true;
}

static void snap_and_satisfy_quality(Mesh* mesh, AdaptOpts const& opts) {
#ifdef OMEGA_H_USE_EGADS
  if (opts.egads_model) {
    ScopedTimer snap_timer("snap");
    mesh->set_parting(OMEGA_H_GHOSTED);
    auto warp = egads_get_snap_warp(
        mesh, opts.egads_model, opts.verbosity >= EACH_REBUILD);
    if (opts.should_smooth_snap) {
      if (opts.verbosity >= EACH_REBUILD) {
        std::cout << "Solving Laplacian of warp field...\n";
      }
      auto t0 = now();
      warp =
          solve_laplacian(mesh, warp, mesh->dim(), opts.snap_smooth_tolerance);
      auto t1 = now();
      if (opts.verbosity >= EACH_REBUILD) {
        std::cout << "Solving Laplacian of warp field took " << (t1 - t0)
                  << " seconds\n";
      }
    }
    mesh->add_tag(VERT, "warp", mesh->dim(), warp);
    while (warp_to_limit(mesh, opts, opts.allow_snap_failure)) {
      if (!satisfy_quality(mesh, opts)) {
        mesh->remove_tag(VERT, "warp");
        break;
      }
    }
  } else
#endif
    satisfy_quality(mesh, opts);
}

static void post_adapt(
    Mesh* mesh, AdaptOpts const& opts, Now t0, Now t1, Now t2, Now t3, Now t4) {
  if (opts.verbosity == EACH_ADAPT) {
    if (!mesh->comm()->rank()) std::cout << "after adapting:\n";
    print_adapt_status(mesh, opts);
  }
  if (opts.verbosity >= EXTRA_STATS) print_adapt_histograms(mesh, opts);
  if (opts.verbosity > SILENT && !mesh->comm()->rank()) {
    std::cout << "addressing edge lengths took " << (t2 - t1) << " seconds\n";
  }
  if (opts.verbosity > SILENT && !mesh->comm()->rank()) {
#ifdef OMEGA_H_USE_EGADS
    if (opts.egads_model) std::cout << "snapping while ";
#endif
    std::cout << "addressing element qualities took " << (t3 - t2);
    std::cout << " seconds\n";
  }
  if (opts.verbosity > SILENT && should_conserve_any(mesh, opts.xfer_opts) &&
      !mesh->comm()->rank()) {
    std::cout << "correcting integral errors took " << (t4 - t3)
              << " seconds\n";
  }
  Now t5 = now();
  if (opts.verbosity > SILENT && !mesh->comm()->rank()) {
    std::cout << "adapting took " << (t5 - t0) << " seconds\n\n";
  }
}

bool adapt(Mesh* mesh, AdaptOpts const& opts) {
  ScopedTimer adapt_timer("adapt");
  OMEGA_H_CHECK(mesh->family() == OMEGA_H_SIMPLEX);
  auto t0 = now();
  if (!pre_adapt(mesh, opts)) return false;
  setup_conservation_tags(mesh, opts);
  auto t1 = now();
  satisfy_lengths(mesh, opts);
  auto t2 = now();
  snap_and_satisfy_quality(mesh, opts);
  auto t3 = now();
  correct_integral_errors(mesh, opts);
  auto t4 = now();
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  post_adapt(mesh, opts, t0, t1, t2, t3, t4);
  return true;
}

}  // end namespace Omega_h
