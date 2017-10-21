#include "Omega_h_adapt.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_shape.hpp"

#include <iostream>

namespace Omega_h {

static void check_okay(Mesh* mesh, AdaptOpts const& opts) {
  auto minq = min_fixable_quality(mesh, opts);
  if (minq < opts.min_quality_allowed) {
    Omega_h_fail("minimum element quality %f < minimum allowed quality %f\n",
        minq, opts.min_quality_allowed);
  }
  auto maxl = mesh->max_length();
  if (maxl > opts.max_length_allowed) {
    Omega_h_fail("maximum edge length %f > maximum allowed length %f\n", maxl,
        opts.max_length_allowed);
  }
}

static bool okay(Mesh* mesh, AdaptOpts const& opts) {
  auto minq = min_fixable_quality(mesh, opts);
  auto maxl = mesh->max_length();
  return minq >= opts.min_quality_allowed && maxl <= opts.max_length_allowed;
}

bool warp_to_limit(Mesh* mesh, AdaptOpts const& opts,
    bool exit_on_stall, Int max_niters) {
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
  do {
    ++i;
    if (i > max_niters) {
      if (exit_on_stall) {
        if (can_print(mesh)) {
          std::cout << "warp_to_limit stalled, dropping warp field and continuing anyway\n";
        }
        mesh->remove_tag(VERT, "warp");
        return true;
      }
      Omega_h_fail(
          "warp step %d : Omega_h is probably unable to satisfy"
          " this warp under this size field\n"
          "min quality %.2e max length %.2e\n",
          i, min_fixable_quality(mesh, opts), mesh->max_length());
    }
    auto half_warp = multiply_each_by(warp, 1.0 / 2.0);
    warp = half_warp;
    remainder = add_each(remainder, half_warp);
    mesh->set_coords(add_each(coords, warp));
  } while (!okay(mesh, opts));
  mesh->set_tag(VERT, "warp", remainder);
  return true;
}

bool approach_metric(Mesh* mesh, AdaptOpts const& opts) {
  auto name = "metric";
  auto target_name = "target_metric";
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
    auto current = interpolate_between_metrics(mesh->nverts(), orig, target, t);
    mesh->set_tag(VERT, name, current);
  } while (!okay(mesh, opts));
  return true;
}

}  // end namespace Omega_h
