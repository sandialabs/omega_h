#include <Omega_h_adapt.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif

#include <cstdio>
#include <iostream>

namespace Omega_h {

template <Int dim>
void compute_ill_metric_dim(
    Mesh* mesh, AdaptOpts const& opts, Omega_h_Isotropy isotropy) {
  mesh->remove_tag(VERT, "metric");
  std::cerr << "getting element implied metrics\n";
  auto elem_metrics = get_element_implied_length_metrics(mesh);
  elem_metrics = apply_isotropy(mesh->nelems(), elem_metrics, isotropy);
  auto ncomps = divide_no_remainder(elem_metrics.size(), mesh->nelems());
  mesh->add_tag(dim, "metric", ncomps, elem_metrics);
  std::cerr << "projecting them\n";
  auto metrics = project_metrics(mesh, elem_metrics);
  OMEGA_H_CHECK(ncomps == divide_no_remainder(metrics.size(), mesh->nverts()));
  if (isotropy == OMEGA_H_ANISOTROPIC) {
    for (Int i = 0; i < 5; ++i) {
      std::cerr << "metric smoothing iteration " << i << '\n';
      metrics = smooth_metric_once(mesh, metrics, true);
    }
  }
  mesh->add_tag(VERT, "metric", ncomps, metrics);
  std::cerr << "computing qualities\n";
  mesh->ask_qualities();
  std::cerr << "writing predicted.vtu\n";
  vtk::write_vtu("predicted.vtu", mesh);
  auto prelim_quals = mesh->get_array<Real>(mesh->dim(), "quality");
  auto min_qual = opts.min_quality_allowed;
  std::cerr << "removing low-quality contributions\n";
  auto elem_metrics_w = deep_copy(elem_metrics);
  auto f = OMEGA_H_LAMBDA(LO e) {
    if (prelim_quals[e] < min_qual) {
      for (Int i = 0; i < ncomps; ++i) {
        elem_metrics_w[e * ncomps + i] = 0.;
      }
    }
  };
  parallel_for(mesh->nelems(), f, "ignore_low_quality_implied");
  elem_metrics = Reals(elem_metrics_w);
  std::cerr << "done removing low-quality contributions\n";
  mesh->remove_tag(dim, "metric");
  mesh->add_tag(dim, "metric", ncomps, elem_metrics);
  std::cerr << "projecting metrics elem -> node\n";
  metrics = project_metrics(mesh, elem_metrics);
  mesh->remove_tag(VERT, "metric");
  mesh->add_tag(VERT, "metric",
      divide_no_remainder(metrics.size(), mesh->nverts()), metrics);
  vtk::write_vtu("no_low_qual_metric.vtu", mesh);
  mesh->remove_tag(VERT, "metric");
  mesh->add_tag(VERT, "metric",
      divide_no_remainder(metrics.size(), mesh->nverts()), metrics);
  vtk::write_vtu("unlimited_ill_metric.vtu", mesh);
  std::cerr << "limiting metric gradation\n";
  metrics = limit_metric_gradation(mesh, metrics, 1.0);
  mesh->remove_tag(VERT, "metric");
  mesh->add_tag(VERT, "metric",
      divide_no_remainder(metrics.size(), mesh->nverts()), metrics);
  mesh->ask_qualities();
  std::cerr << "writing corrected.vtu\n";
  vtk::write_vtu("corrected.vtu", mesh);
  std::cerr << "done with \"ill\" metric\n";
}

static void compute_ill_metric(
    Mesh* mesh, AdaptOpts const& opts, Omega_h_Isotropy isotropy) {
  if (mesh->dim() == 3)
    compute_ill_metric_dim<3>(mesh, opts, isotropy);
  else if (mesh->dim() == 2)
    compute_ill_metric_dim<2>(mesh, opts, isotropy);
  else if (mesh->dim() == 1)
    compute_ill_metric_dim<1>(mesh, opts, isotropy);
  else
    OMEGA_H_NORETURN();
}

static void compute_typical_metric(Mesh* mesh, Omega_h_Isotropy isotropy) {
  auto metrics = get_implied_metrics(mesh);
  metrics = apply_isotropy(mesh->nverts(), metrics, isotropy);
  metrics = limit_metric_gradation(mesh, metrics, 1.0);
  auto ncomps = divide_no_remainder(metrics.size(), mesh->nverts());
  mesh->remove_tag(VERT, "metric");
  mesh->add_tag(VERT, "metric", ncomps, metrics);
}

void fix_for_given_metric(
    Mesh* mesh, AdaptOpts const& adapt_opts, bool verbose) {
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
              << "maxlen(" << maxlen << "), max_length_desired*2("
              << opts.max_length_desired * 2.0 << "))\n";
    opts.verbosity = EXTRA_STATS;
    opts.nsliver_layers = 10;
    opts.min_quality_desired = min2(minqual + 0.1, 1.0);
#ifdef OMEGA_H_USE_EGADS
    opts.allow_snap_failure = true;
#endif
    adapt(mesh, opts);
    minqual = mesh->min_quality();
    maxlen = mesh->max_length();
    if (verbose) std::cout << "minimum quality " << minqual << '\n';
    if (verbose) std::cout << "maximum length " << maxlen << '\n';
    if (minqual < minqual_old + 1e-3) break;  // stalled
  }
}

void fix(Mesh* mesh, AdaptOpts const& adapt_opts, Omega_h_Isotropy isotropy,
    bool verbose) {
  verbose = verbose && can_print(mesh);
  if (verbose) std::cout << "Computing bad-mesh \"implied\" metric\n";
  compute_ill_metric(mesh, adapt_opts, isotropy);
  fix_for_given_metric(mesh, adapt_opts, verbose);
  // TODO: commenting this out because metric quality corrections assume
  // isotropy!!! we need to fix this, but for now just skip this step
  if ((0)) {
    if (verbose) std::cout << "Computing more typical \"implied\" metric\n";
    compute_typical_metric(mesh, isotropy);
    fix_for_given_metric(mesh, adapt_opts, verbose);
  }
}

void grade_fix_adapt(
    Mesh* mesh, AdaptOpts const& opts, Reals target_metric, bool verbose) {
  verbose = verbose && can_print(mesh);
  auto metric_dim = get_metrics_dim(mesh->nverts(), target_metric);
  if (verbose) std::cout << "Limiting target metric gradation...\n";
  target_metric = Omega_h::limit_metric_gradation(mesh, target_metric, 0.5);
  mesh->remove_tag(VERT, "target_metric");
  mesh->add_tag(VERT, "target_metric", symm_ncomps(metric_dim), target_metric);
  if (mesh->has_tag(0, "metric")) {
    if (verbose)
      std::cout << "Mesh already has \"metric\" on vertices, assuming that is "
                   "current\n";
  } else {
    if (verbose) std::cout << "Deriving implied metric...\n";
    if (metric_dim == 1) {
      Omega_h::add_implied_isos_tag(mesh);
    } else {
      Omega_h::add_implied_metric_tag(mesh);
    }
  }
  auto min_qual = mesh->min_quality();
  if (verbose) std::cout << "Initial mesh has minimum quality " << min_qual;
  if (min_qual < opts.min_quality_allowed) {
    if (verbose) {
      std::cout << " < minimum acceptable quality " << opts.min_quality_allowed
                << '\n';
      std::cout
          << "Omega_h will now attempt to repair the initial mesh quality.\n";
      std::cout << "This could take some time...\n";
    }
    auto isotropy =
        (metric_dim == 1 ? OMEGA_H_ISO_LENGTH : OMEGA_H_ANISOTROPIC);
    Omega_h::fix(mesh, opts, isotropy, /*verbose=*/true);
    if (verbose) std::cout << "\nOmega_h is done repairing mesh quality!\n\n";
  } else {
    if (verbose) std::cout << ", which is good\n";
  }
  if (verbose) std::cout << "Adapting...\n";
  while (Omega_h::approach_metric(mesh, opts)) {
    Omega_h::adapt(mesh, opts);
  }
  if (verbose) std::cout << "\nDone adapting!\n\n";
}

}  // namespace Omega_h
