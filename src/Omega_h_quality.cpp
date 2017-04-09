#include "Omega_h_quality.hpp"

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
Reals measure_qualities_tmpl(Mesh* mesh, LOs a2e) {
  MetricElementQualities<mesh_dim, metric_dim> measurer(mesh);
  auto ev2v = mesh->ask_verts_of(mesh_dim);
  auto na = a2e.size();
  Write<Real> qualities(na);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<mesh_dim + 1>(ev2v, e);
    qualities[a] = measurer.measure(v);
  };
  parallel_for(na, f);
  return qualities;
}

Reals measure_qualities(Mesh* mesh, LOs a2e) {
  auto metrics = mesh->get_array<Real>(VERT, "metric");
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return measure_qualities_tmpl<3, 3>(mesh, a2e);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return measure_qualities_tmpl<2, 2>(mesh, a2e);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return measure_qualities_tmpl<3, 1>(mesh, a2e);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return measure_qualities_tmpl<2, 1>(mesh, a2e);
  }
  if (mesh->dim() == 1) {
    return Reals(a2e.size(), 1.0);
  }
  NORETURN(Reals());
}

Reals measure_qualities(Mesh* mesh) {
  return measure_qualities(mesh, LOs(mesh->nelems(), 0, 1));
}

}  // end namespace Omega_h
