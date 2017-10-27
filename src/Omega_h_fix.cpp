#include <Omega_h_adapt.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_file.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif

#include <iostream>
#include <cstdio>

namespace Omega_h {

static void compute_ill_metric(Mesh* mesh) {
  auto elem_metrics = get_element_implied_length_metrics(mesh);
  auto metrics = project_metrics(mesh, elem_metrics);
  metrics = limit_metric_gradation(mesh, metrics, 1.0);
  mesh->remove_tag(VERT, "metric");
  mesh->add_tag(
      VERT, "metric", symm_ncomps(mesh->dim()), metrics);
}

static void compute_typical_metric(Mesh* mesh, Omega_h_Isotropy isotropy) {
  auto metrics = get_implied_metrics(mesh);
  metrics = apply_isotropy(mesh->nverts(), metrics, isotropy);
  auto ncomps = divide_no_remainder(metrics.size(), mesh->nverts());
  mesh->remove_tag(VERT, "metric");
  mesh->add_tag(VERT, "metric", ncomps, metrics);
}

static void fix_for_given_metric(Mesh* mesh, AdaptOpts const& adapt_opts, bool verbose) {
  if (verbose) std::cout << "computing minimum quality\n";
  auto minqual = mesh->min_quality();
  auto maxlen = mesh->max_length();
  if (verbose) std::cout << "minimum quality " << minqual << '\n';
  if (verbose) std::cout << "maximum length " << maxlen << '\n';
  AdaptOpts opts = adapt_opts;
  while (true) {
    auto minqual_old = minqual;
    opts.min_quality_allowed = minqual;
    opts.max_length_allowed = max2(maxlen, opts.max_length_desired * 2.0);
    std::cout << "max_length_allowed(" << opts.max_length_allowed << ") = max("
      << "maxlen(" << maxlen << "), max_length_desired*2(" << opts.max_length_desired * 2.0 << "))\n";
    opts.verbosity = EXTRA_STATS;
    opts.nsliver_layers = 10;
    opts.min_quality_desired = min2(minqual + 0.1, 1.0);
#ifdef OMEGA_H_USE_EGADS
    opts.allow_snap_failure = true;
#endif
    adapt(mesh, opts);
    if ((0)) {
      std::cout << "writing debug.osh after adapt\n";
      Omega_h::binary::write("debug.osh", mesh);
    }
    minqual = mesh->min_quality();
    maxlen = mesh->max_length();
    if (verbose) std::cout << "minimum quality " << minqual << '\n';
    if (verbose) std::cout << "maximum length " << maxlen << '\n';
    if (minqual < minqual_old + 1e-3) break;  // stalled
  }
}

void fix(Mesh* mesh
    , AdaptOpts const& adapt_opts
    , Omega_h_Isotropy isotropy
    , bool verbose
    ) {
  verbose = verbose && can_print(mesh);
  if (verbose) std::cout << "Computing bad-mesh \"implied\" metric\n";
  compute_ill_metric(mesh);
  fix_for_given_metric(mesh, adapt_opts, verbose);
  if (verbose) std::cout << "Computing more typical \"implied\" metric\n";
  compute_typical_metric(mesh, isotropy);
  fix_for_given_metric(mesh, adapt_opts, verbose);
}

}
