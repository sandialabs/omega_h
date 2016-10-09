#include "adapt.hpp"

#include <iomanip>
#include <iostream>

#include "array.hpp"
#include "coarsen.hpp"
#include "mark.hpp"
#include "quality.hpp"
#include "refine.hpp"
#include "simplices.hpp"
#include "swap.hpp"
#include "timer.hpp"

namespace Omega_h {

static void goal_stats(Mesh* mesh, char const* name, Int ent_dim, Reals values,
    Real floor, Real ceil, Real minval, Real maxval) {
  auto low_marks = each_lt(values, floor);
  auto high_marks = each_gt(values, ceil);
  auto nlow = count_owned_marks(mesh, ent_dim, low_marks);
  auto nhigh = count_owned_marks(mesh, ent_dim, high_marks);
  auto ntotal = mesh->nglobal_ents(ent_dim);
  auto nmid = ntotal - nlow - nhigh;
  if (mesh->comm()->rank() == 0) {
    auto precision_before = std::cout.precision();
    std::ios::fmtflags stream_state(std::cout.flags());
    std::cout << std::fixed << std::setprecision(2);
    std::cout << ntotal << " " << plural_names[ent_dim];
    std::cout << ", " << name << " [" << minval << "," << maxval << "]";
    if (nlow) {
      std::cout << ", " << nlow << " <" << floor;
    }
    if (nmid) {
      std::cout << ", " << nmid << " in [" << floor << "," << ceil << "]";
    }
    if (nhigh) {
      std::cout << ", " << nhigh << " >" << ceil;
    }
    std::cout << '\n';
    std::cout.flags(stream_state);
    std::cout.precision(precision_before);
  }
}

static void get_minmax(
    Mesh* mesh, Reals values, Real* p_minval, Real* p_maxval) {
  *p_minval = mesh->comm()->allreduce(min(values), OMEGA_H_MIN);
  *p_maxval = mesh->comm()->allreduce(max(values), OMEGA_H_MAX);
}

static void adapt_summary(Mesh* mesh, Real qual_floor, Real qual_ceil,
    Real minqual, Real maxqual, Real minlen,
    Real maxlen) {
  goal_stats(mesh, "quality", mesh->dim(), mesh->ask_qualities(), qual_floor,
      qual_ceil, minqual, maxqual);
  goal_stats(mesh, "length", EDGE, mesh->ask_lengths(),
      min_length_desired, max_length_desired, minlen, maxlen);
}

bool adapt_check(Mesh* mesh, AdaptOpts const& opts) {
  Real minqual, maxqual;
  get_minmax(mesh, mesh->ask_qualities(), &minqual, &maxqual);
  Real minlen, maxlen;
  get_minmax(mesh, mesh->ask_lengths(), &minlen, &maxlen);
  if (minqual >= opts.minimum_quality_desired &&
      minlen >= min_length_desired &&
      maxlen <= max_length_desired) {
    if (opts.verbosity > SILENT && mesh->comm()->rank() == 0) {
      std::cout << "mesh is good: quality [" << minqual << "," << maxqual
                << "], length [" << minlen << "," << maxlen << "]\n";
    }
    return true;
  }
  if (opts.verbosity > SILENT) {
    adapt_summary(mesh, qual_floor, qual_ceil, minqual,
        maxqual, minlen, maxlen);
  }
  return false;
}

static void do_histogram(Mesh* mesh) {
  auto histogram = get_quality_histogram(mesh);
  if (mesh->comm()->rank() == 0) print_quality_histogram(histogram);
}

static void validate(Mesh* mesh, AdaptOpts const& opts) {
  CHECK(0.0 <= opts.minimum_quality_allowed);
  CHECK(opts.minimum_quality_allowed <=
        opts.minimum_quality_desired);
  CHECK(opts.minimum_quality_desired <= 1.0);
  CHECK(opts.nsliver_layers_considered >= 0);
  CHECK(opts.nsliver_layers_considered < 100);
  CHECK(mesh->min_quality() >= opts.minimum_quality_allowed);
}

static bool pre_adapt(Mesh* mesh, AdaptOpts const& opts) {
  validate(mesh, opts);
  if (opts.verbosity >= EACH_ADAPT && !mesh->comm->rank()) {
    std::cout << "before adapting:\n";
  }
  if (adapt_check(mesh, opts)) return false;
  if (verbosity >= EXTRA_STATS) do_histogram(mesh);
  if ((verbosity >= EACH_REBUILD) && mesh->comm->rank() == 0) {
    std::cout << "addressing edge lengths\n";
  }
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
  if ((verbosity >= EACH_REBUILD) && !mesh->comm->rank()) {
    std::cout << "addressing element qualities\n";
  }
  do {
    if (swap_edges(mesh, qual_ceil, nlayers, (verbosity >= 2))) {
      post_rebuild(mesh, opts);
      continue;
    }
    if (coarsen_slivers(mesh, opts)) {
      post_rebuild(mesh, opts);
      continue;
    }
    if ((verbosity > SILENT) && mesh->comm->rank() == 0) {
      std::cout << "adapt() could not satisfy quality\n";
    }
    break;
  } while (mesh->min_quality() < opts.min_quality_desired);
}

static void post_adapt(Mesh* mesh, AdaptOpts const& opts,
    Now t0, Now t1, Now t2, Now t3) {
  if (opts.verbosity == EACH_ADAPT) {
    if (mesh->comm->rank() == 0) std::cout << "after adapting:\n";
    adapt_check(mesh, opts);
  }
  if (verbosity >= EXTRA_STATS) do_histogram(mesh);
  if (verbosity > SILENT && mesh->comm->rank() == 0) {
    std::cout << "addressing edge lengths took " << (t2 - t1) << " seconds\n";
  }
  if (verbosity > SILENT && mesh->comm->rank() == 0) {
    std::cout << "addressing element qualities took " << (t3 - t2)
              << " seconds\n";
  }
  Now t4 = now();
  if (verbosity > SILENT && mesh->comm->rank() == 0) {
    std::cout << "adapting took " << (t4 - t0) << " seconds\n\n";
  }
}

bool adapt(Mesh* mesh, AdaptOpts const& opts) {
  Now t0 = now();
  pre_adapt(mesh, opts);
  Now t1 = now();
  satisfy_lengths(mesh, opts);
  Now t2 = now();
  satisfy_quality(mesh, opts);
  Now t3 = now();
  post_adapt(mesh, opts, t0, t1, t2, t3);
  return true;
}

}  // end namespace Omega_h
