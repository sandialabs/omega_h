#include <iomanip>
#include <iostream>

#include "array.hpp"
#include "coarsen.hpp"
#include "control.hpp"
#include "histogram.hpp"
#include "quality.hpp"
#include "refine.hpp"
#include "swap.hpp"
#include "timer.hpp"
#include "laplace.hpp"

#ifdef OMEGA_H_USE_EGADS
#include "Omega_h_egads.hpp"
#endif

namespace Omega_h {

AdaptOpts::AdaptOpts(Mesh* mesh) {
  min_length_desired = 1.0 / sqrt(2.0);
  max_length_desired = sqrt(2.0);
  max_length_allowed = ArithTraits<Real>::max();
  if (mesh->dim() == 3) {
    min_quality_allowed = 0.20;
    min_quality_desired = 0.30;
  }
  if (mesh->dim() == 2) {
    min_quality_allowed = 0.30;
    min_quality_desired = 0.40;
  }
  nsliver_layers = 4;
  verbosity = EACH_REBUILD;
  length_histogram_min = 0.0;
  length_histogram_max = 3.0;
#ifdef OMEGA_H_USE_EGADS
  egads_model = nullptr;
  should_smooth_snap = true;
  snap_smooth_tolerance = 1e-2;
#endif
}

static void adapt_summary(Mesh* mesh, AdaptOpts const& opts,
    MinMax<Real> qualstats, MinMax<Real> lenstats) {
  print_goal_stats(mesh, "quality", mesh->dim(), mesh->ask_qualities(),
      {opts.min_quality_allowed, opts.min_quality_desired}, qualstats);
  print_goal_stats(mesh, "length", EDGE, mesh->ask_lengths(),
      {opts.min_length_desired, opts.max_length_desired}, lenstats);
}

static bool adapt_check(Mesh* mesh, AdaptOpts const& opts) {
  auto qualstats = get_minmax(mesh->comm(), mesh->ask_qualities());
  auto lenstats = get_minmax(mesh->comm(), mesh->ask_lengths());
  if (qualstats.min >= opts.min_quality_desired &&
      lenstats.min >= opts.min_length_desired &&
      lenstats.max <= opts.max_length_desired) {
    if (opts.verbosity > SILENT && mesh->comm()->rank() == 0) {
      std::cout << "mesh is good: quality [" << qualstats.min << "," << qualstats.max
                << "], length [" << lenstats.min << "," << lenstats.max << "]\n";
    }
    return true;
  }
  if (opts.verbosity > SILENT) {
    adapt_summary(mesh, opts, qualstats, lenstats);
  }
  return false;
}

static void do_histograms(Mesh* mesh, AdaptOpts const& opts) {
  auto qh =
      get_histogram<10>(mesh, mesh->dim(), mesh->ask_qualities(), 0.0, 1.0);
  print_histogram(mesh, qh, "quality");
  auto lh = get_histogram<10>(mesh, VERT, mesh->ask_lengths(),
      opts.length_histogram_min, opts.length_histogram_max);
  print_histogram(mesh, lh, "length");
}

static void validate(Mesh* mesh, AdaptOpts const& opts) {
  CHECK(0.0 <= opts.min_quality_allowed);
  CHECK(opts.min_quality_allowed <= opts.min_quality_desired);
  CHECK(opts.min_quality_desired <= 1.0);
  CHECK(opts.nsliver_layers >= 0);
  CHECK(opts.nsliver_layers < 100);
  auto mq = mesh->min_quality();
  if (mq < opts.min_quality_allowed && !mesh->comm()->rank()) {
    std::cout << "WARNING: worst input element has quality " << mq
              << " but minimum allowed is " << opts.min_quality_allowed << "\n";
  }
}

static bool pre_adapt(Mesh* mesh, AdaptOpts const& opts) {
  validate(mesh, opts);
  if (opts.verbosity >= EACH_ADAPT && !mesh->comm()->rank()) {
    std::cout << "before adapting:\n";
  }
  if (adapt_check(mesh, opts)) return false;
  if (opts.verbosity >= EXTRA_STATS) do_histograms(mesh, opts);
  if ((opts.verbosity >= EACH_REBUILD) && !mesh->comm()->rank()) {
    std::cout << "addressing edge lengths\n";
  }
  return true;
}

static void post_rebuild(Mesh* mesh, AdaptOpts const& opts) {
  if (opts.verbosity >= EACH_REBUILD) adapt_check(mesh, opts);
}

static void satisfy_lengths(Mesh* mesh, AdaptOpts const& opts) {
  bool did_anything;
  do {
    did_anything = false;
    if (refine_by_size(mesh, opts)) {
      post_rebuild(mesh, opts);
      did_anything = true;
    }
    if (coarsen_by_size(mesh, opts)) {
      post_rebuild(mesh, opts);
      did_anything = true;
    }
  } while (did_anything);
}

static void satisfy_quality(Mesh* mesh, AdaptOpts const& opts) {
  if (mesh->min_quality() >= opts.min_quality_desired) return;
  if ((opts.verbosity >= EACH_REBUILD) && !mesh->comm()->rank()) {
    std::cout << "addressing element qualities\n";
  }
  do {
    if (swap_edges(mesh, opts)) {
      post_rebuild(mesh, opts);
      continue;
    }
    if (coarsen_slivers(mesh, opts)) {
      post_rebuild(mesh, opts);
      continue;
    }
    if ((opts.verbosity > SILENT) && !mesh->comm()->rank()) {
      std::cout << "adapt() could not satisfy quality\n";
    }
    break;
  } while (mesh->min_quality() < opts.min_quality_desired);
}

static void snap_and_satisfy_quality(Mesh* mesh, AdaptOpts const& opts) {
#ifdef OMEGA_H_USE_EGADS
  if (opts.egads_model) {
    mesh->set_parting(OMEGA_H_GHOSTED);
    auto warp = egads_get_snap_warp(mesh, opts.egads_model);
    if (opts.should_smooth_snap) {
      warp = solve_laplacian(mesh, warp, mesh->dim(), opts.snap_smooth_tolerance);
    }
    mesh->add_tag(VERT, "warp", mesh->dim(), OMEGA_H_LINEAR_INTERP,
        OMEGA_H_DO_OUTPUT, warp);
    while (warp_to_limit(mesh, opts)) satisfy_quality(mesh, opts);
  } else
#endif
  satisfy_quality(mesh, opts);
}

static void post_adapt(
    Mesh* mesh, AdaptOpts const& opts, Now t0, Now t1, Now t2, Now t3) {
  if (opts.verbosity == EACH_ADAPT) {
    if (!mesh->comm()->rank()) std::cout << "after adapting:\n";
    adapt_check(mesh, opts);
  }
  if (opts.verbosity >= EXTRA_STATS) do_histograms(mesh, opts);
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
  Now t4 = now();
  if (opts.verbosity > SILENT && !mesh->comm()->rank()) {
    std::cout << "adapting took " << (t4 - t0) << " seconds\n\n";
    add_to_global_timer("adapting", t4 - t0);
  }
}

bool adapt(Mesh* mesh, AdaptOpts const& opts) {
  auto t0 = now();
  if (!pre_adapt(mesh, opts)) return false;
  auto t1 = now();
  satisfy_lengths(mesh, opts);
  auto t2 = now();
  snap_and_satisfy_quality(mesh, opts);
  auto t3 = now();
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  post_adapt(mesh, opts, t0, t1, t2, t3);
  return true;
}

}  // end namespace Omega_h
