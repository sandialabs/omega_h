#include "array.hpp"
#include "internal.hpp"
#include "metric.hpp"
#include "size.hpp"

namespace Omega_h {

bool warp_to_limit(Mesh* mesh, AdaptOpts const& opts) {
  auto min_qual = opts.min_quality_allowed;
  if (!mesh->has_tag(VERT, "warp")) return false;
  CHECK(mesh->min_quality() >= min_qual);
  auto coords = mesh->coords();
  auto warp = mesh->get_array<Real>(VERT, "warp");
  mesh->set_coords(add_each(coords, warp));
  if (mesh->min_quality() >= min_qual) {
    mesh->remove_tag(VERT, "warp");
    return true;
  }
  auto remainder = Reals(warp.size(), 0.0);
  do {
    auto half_warp = multiply_each_by(1.0 / 2.0, warp);
    warp = half_warp;
    remainder = add_each(remainder, half_warp);
    mesh->set_coords(add_each(coords, warp));
  } while (mesh->min_quality() < min_qual);
  mesh->set_tag(VERT, "warp", remainder);
  return true;
}

static bool okay(Mesh* mesh, AdaptOpts const& opts) {
  return mesh->min_quality() >= opts.min_quality_allowed &&
         mesh->max_length() <= opts.max_length_desired * 2.0;
}

static bool approach_either(Mesh* mesh, AdaptOpts const& opts, std::string const& name,
    Reals (*interpolator)(Int dim, Reals orig, Reals target, Real t)) {
  auto target_name = std::string("target_") + name;
  if (!mesh->has_tag(VERT, target_name)) return false;
  CHECK(okay(mesh, opts));
  auto orig = mesh->get_array<Real>(VERT, name);
  auto target = mesh->get_array<Real>(VERT, target_name);
  mesh->set_tag(VERT, name, target);
  if (okay(mesh, opts)) {
    mesh->remove_tag(VERT, target_name);
    return true;
  }
  Real t = 1.0;
  do {
    t /= 2.0;
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
    return approach_either(
        mesh, opts, "metric", &interpolate_between_metrics);
  }
  NORETURN(true);
}

}  // end namespace Omega_h
