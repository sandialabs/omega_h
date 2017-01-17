#include "array.hpp"
#include "internal.hpp"
#include "metric.hpp"
#include "size.hpp"

#include <cstdio>
#include <iostream>

namespace Omega_h {

static void check_okay(Mesh* mesh, AdaptOpts const& opts) {
  auto minq = mesh->min_quality();
  if (minq < opts.min_quality_allowed) {
    Omega_h_fail("minimum element quality %f < minimum allowed quality %f\n",
        minq, opts.min_quality_allowed);
  }
  auto maxl = mesh->max_length();
  if (maxl > opts.max_length_allowed) {
    Omega_h_fail("maximum edge length %f > maximum allowed length %f\n", maxl,
        opts.max_length_allowed);
  }
  fprintf(stderr, "check_okay() minqual %.2e maxlen %.2e\n", minq, maxl);
}

static bool okay(Mesh* mesh, AdaptOpts const& opts) {
  auto minq = mesh->min_quality();
  auto maxl = mesh->max_length();
  fprintf(stderr, "okay() minqual %.2e maxlen %.2e\n", minq, maxl);
  return minq >= opts.min_quality_allowed && maxl <= opts.max_length_allowed;
}

bool warp_to_limit(Mesh* mesh, AdaptOpts const& opts) {
  std::cerr << "warp_to_limit()\n";
  if (!mesh->has_tag(VERT, "warp")) return false;
  check_okay(mesh, opts);
  auto coords = mesh->coords();
  auto warp = mesh->get_array<Real>(VERT, "warp");
  mesh->set_coords(add_each(coords, warp));
  if (okay(mesh, opts)) {
    mesh->remove_tag(VERT, "warp");
    return true;
  }
  auto remainder = Reals(warp.size(), 0.0);
  Int i = 0;
  constexpr Int max_i = 40;
  do {
    ++i;
    if (i == max_i) {
      vtk::write_vtu("warp_fail.vtu", mesh, mesh->dim());
      vtk::write_vtu("warp_fail_edges.vtu", mesh, EDGE);
      Omega_h_fail("warp step %d : Omega_h is probably unable to satisfy"
                   " this warp under this size field\n"
                   "min quality %.2e max length %.2e\n",
                   i, mesh->min_quality(), mesh->max_length());
    }
    auto half_warp = multiply_each_by(1.0 / 2.0, warp);
    warp = half_warp;
    remainder = add_each(remainder, half_warp);
    mesh->set_coords(add_each(coords, warp));
  } while (!okay(mesh, opts));
  mesh->set_tag(VERT, "warp", remainder);
  return true;
}

static bool approach_either(Mesh* mesh, AdaptOpts const& opts,
    std::string const& name,
    Reals (*interpolator)(Int dim, Reals orig, Reals target, Real t)) {
  auto target_name = std::string("target_") + name;
  if (!mesh->has_tag(VERT, target_name)) return false;
  check_okay(mesh, opts);
  auto orig = mesh->get_array<Real>(VERT, name);
  auto target = mesh->get_array<Real>(VERT, target_name);
  mesh->set_tag(VERT, name, target);
  if (okay(mesh, opts)) {
    mesh->remove_tag(VERT, target_name);
    return true;
  }
  Real t = 1.0;
  constexpr Real min_t = 1e-4;
  do {
    t /= 2.0;
    if (t < min_t) {
      Omega_h_fail(
          "size field approach step = %f < %f.\n"
          "Omega_h is probably unable to satisfy this size field\n",
          t, min_t);
    }
    auto current = (*interpolator)(mesh->dim(), orig, target, t);
    mesh->set_tag(VERT, name, current);
  } while (!okay(mesh, opts));
  return true;
}

static Reals isos_wrapper(Int, Reals orig, Reals target, Real t) {
  return interpolate_between_isos(orig, target, t);
}

bool approach_size_field(Mesh* mesh, AdaptOpts const& opts) {
  if (mesh->has_tag(VERT, "size")) {
    return approach_either(mesh, opts, "size", &isos_wrapper);
  }
  if (mesh->has_tag(VERT, "metric")) {
    return approach_either(mesh, opts, "metric", &interpolate_between_metrics);
  }
  NORETURN(true);
}

}  // end namespace Omega_h
