#include "quality.hpp"

#include <iomanip>
#include <iostream>

#include "array.hpp"
#include "loop.hpp"
#include "map.hpp"

namespace Omega_h {

template <typename ElementQualities, Int dim>
Reals measure_qualities_tmpl(Mesh* mesh, LOs a2e) {
  ElementQualities measurer(mesh);
  auto ev2v = mesh->ask_verts_of(mesh->dim());
  auto na = a2e.size();
  Write<Real> qualities(na);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<dim + 1>(ev2v, e);
    qualities[a] = measurer.measure(v);
  };
  parallel_for(na, f);
  return qualities;
}

Reals measure_qualities(Mesh* mesh, LOs a2e) {
  if (mesh->dim() == 3) {
    if (mesh->has_tag(VERT, "metric")) {
      return measure_qualities_tmpl<MetricElementQualities, 3>(mesh, a2e);
    } else {
      return measure_qualities_tmpl<RealElementQualities, 3>(mesh, a2e);
    }
  } else {
    CHECK(mesh->dim() == 2);
    if (mesh->has_tag(VERT, "metric")) {
      return measure_qualities_tmpl<MetricElementQualities, 2>(mesh, a2e);
    } else {
      return measure_qualities_tmpl<RealElementQualities, 2>(mesh, a2e);
    }
  }
}

Reals measure_qualities(Mesh* mesh) {
  return measure_qualities(mesh, LOs(mesh->nelems(), 0, 1));
}

}  // end namespace Omega_h
