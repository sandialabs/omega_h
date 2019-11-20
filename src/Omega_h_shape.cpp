#include "Omega_h_shape.hpp"

#include "Omega_h_for.hpp"

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
Reals measure_edges_metric_tmpl(Mesh* mesh, LOs a2e, Reals metrics) {
  MetricEdgeLengths<mesh_dim, metric_dim> measurer(mesh, metrics);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto na = a2e.size();
  Write<Real> lengths(na);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<2>(ev2v, e);
    lengths[a] = measurer.measure(v);
  };
  parallel_for(na, f, "measure_edges");
  return lengths;
}

Reals measure_edges_metric(Mesh* mesh, LOs a2e, Reals metrics) {
  if (a2e.size() == 0) return Reals({});
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return measure_edges_metric_tmpl<3, 3>(mesh, a2e, metrics);
  } else if (mesh->dim() == 2 && metric_dim == 2) {
    return measure_edges_metric_tmpl<2, 2>(mesh, a2e, metrics);
  } else if (mesh->dim() == 3 && metric_dim == 1) {
    return measure_edges_metric_tmpl<3, 1>(mesh, a2e, metrics);
  } else if (mesh->dim() == 2 && metric_dim == 1) {
    return measure_edges_metric_tmpl<2, 1>(mesh, a2e, metrics);
  } else if (mesh->dim() == 1 && metric_dim == 1) {
    return measure_edges_metric_tmpl<1, 1>(mesh, a2e, metrics);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals measure_edges_metric(Mesh* mesh, LOs a2e) {
  return measure_edges_metric(mesh, a2e, mesh->get_array<Real>(VERT, "metric"));
}

Reals measure_edges_metric(Mesh* mesh, Reals metrics) {
  return measure_edges_metric(mesh, LOs(mesh->nedges(), 0, 1), metrics);
}

Reals measure_edges_metric(Mesh* mesh) {
  return measure_edges_metric(mesh, mesh->get_array<Real>(VERT, "metric"));
}

template <Int sdim, Int edim>
Reals measure_ents_real_tmpl(Mesh* mesh, LOs a2e, Reals coords) {
  OMEGA_H_TIME_FUNCTION;
  RealSimplexSizes measurer(coords);
  auto ev2v = mesh->ask_verts_of(edim);
  auto na = a2e.size();
  Write<Real> sizes(na);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<edim + 1>(ev2v, e);
    sizes[a] = measurer.measure<sdim, edim>(v);
  };
  parallel_for(na, f, "measure_ents_real");
  return sizes;
}

Reals measure_ents_real(Mesh* mesh, Int ent_dim, LOs a2e, Reals coords) {
  if (mesh->dim() == 3 && ent_dim == 3)
    return measure_ents_real_tmpl<3, 3>(mesh, a2e, coords);
  if (mesh->dim() == 3 && ent_dim == 2)
    return measure_ents_real_tmpl<3, 2>(mesh, a2e, coords);
  if (mesh->dim() == 3 && ent_dim == 1)
    return measure_ents_real_tmpl<3, 1>(mesh, a2e, coords);
  if (mesh->dim() == 2 && ent_dim == 2)
    return measure_ents_real_tmpl<2, 2>(mesh, a2e, coords);
  if (mesh->dim() == 2 && ent_dim == 1)
    return measure_ents_real_tmpl<2, 1>(mesh, a2e, coords);
  if (mesh->dim() == 1 && ent_dim == 1)
    return measure_ents_real_tmpl<1, 1>(mesh, a2e, coords);
  OMEGA_H_NORETURN(Reals());
}

Reals measure_elements_real(Mesh* mesh, LOs a2e) {
  return measure_ents_real(mesh, mesh->dim(), a2e, mesh->coords());
}

Reals measure_elements_real(Mesh* mesh) {
  return measure_elements_real(mesh, LOs(mesh->nelems(), 0, 1));
}

Reals measure_edges_real(Mesh* mesh, LOs a2e) {
  return measure_ents_real(mesh, EDGE, a2e, mesh->coords());
}

Reals measure_edges_real(Mesh* mesh) {
  return measure_edges_real(mesh, LOs(mesh->nedges(), 0, 1));
}

}  // end namespace Omega_h
