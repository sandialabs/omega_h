#include "array.hpp"
#include "internal.hpp"
#include "metric.hpp"

namespace osh {

bool warp_to_limit(Mesh* mesh, Real min_qual) {
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

bool approach_metric(Mesh* mesh, Real min_qual) {
  if (!mesh->has_tag(VERT, "target_metric")) return false;
  CHECK(mesh->min_quality() >= min_qual);
  auto orig = mesh->get_array<Real>(VERT, "metric");
  auto target = mesh->get_array<Real>(VERT, "target_metric");
  mesh->set_tag(VERT, "metric", target);
  if (mesh->min_quality() >= min_qual) {
    mesh->remove_tag(VERT, "target_metric");
    return true;
  }
  Real t = 1.0;
  do {
    t /= 2.0;
    auto current = interpolate_metrics(mesh->dim(), orig, target, t);
    mesh->set_tag(VERT, "metric", current);
  } while (mesh->min_quality() < min_qual);
  return true;
}

}  // end namespace osh
