#include "Omega_h_shape.hpp"

#include "Omega_h_loop.hpp"

namespace Omega_h {

template <typename EdgeLengths>
static Reals measure_edges_tmpl(Mesh* mesh, LOs a2e, EdgeLengths impl) {
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto na = a2e.size();
  Write<Real> lengths(na);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<2>(ev2v, e);
    lengths[a] = impl.measure(v);
  };
  parallel_for(na, f);
  return lengths;
}

Reals measure_edges_real(Mesh* mesh, LOs a2e) {
  if (mesh->dim() == 3) {
    return measure_edges_tmpl(mesh, a2e, RealEdgeLengths<3>(mesh));
  }
  if (mesh->dim() == 2) {
    return measure_edges_tmpl(mesh, a2e, RealEdgeLengths<2>(mesh));
  }
  if (mesh->dim() == 1) {
    return measure_edges_tmpl(mesh, a2e, RealEdgeLengths<1>(mesh));
  }
  OMEGA_H_NORETURN(Reals());
}

Reals measure_edges_metric(Mesh* mesh, LOs a2e, Reals metrics) {
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return measure_edges_tmpl(
        mesh, a2e, MetricEdgeLengths<3, 3>(mesh, metrics));
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return measure_edges_tmpl(
        mesh, a2e, MetricEdgeLengths<2, 2>(mesh, metrics));
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return measure_edges_tmpl(
        mesh, a2e, MetricEdgeLengths<3, 1>(mesh, metrics));
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return measure_edges_tmpl(
        mesh, a2e, MetricEdgeLengths<2, 1>(mesh, metrics));
  }
  if (mesh->dim() == 1 && metric_dim == 1) {
    return measure_edges_tmpl(
        mesh, a2e, MetricEdgeLengths<1, 1>(mesh, metrics));
  }
  OMEGA_H_NORETURN(Reals());
}

Reals measure_edges_metric(Mesh* mesh, LOs a2e) {
  return measure_edges_metric(mesh, a2e, mesh->get_array<Real>(VERT, "metric"));
}

Reals measure_edges_real(Mesh* mesh) {
  return measure_edges_real(mesh, LOs(mesh->nedges(), 0, 1));
}

Reals measure_edges_metric(Mesh* mesh, Reals metrics) {
  return measure_edges_metric(mesh, LOs(mesh->nedges(), 0, 1), metrics);
}

Reals measure_edges_metric(Mesh* mesh) {
  return measure_edges_metric(mesh, mesh->get_array<Real>(VERT, "metric"));
}

template <Int dim>
static Reals measure_elements_real_tmpl(Mesh* mesh, LOs a2e) {
  RealElementSizes measurer(mesh);
  auto ev2v = mesh->ask_elem_verts();
  auto na = a2e.size();
  Write<Real> sizes(na);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<dim + 1>(ev2v, e);
    sizes[a] = measurer.measure(v);
  };
  parallel_for(na, f);
  return sizes;
}

Reals measure_elements_real(Mesh* mesh, LOs a2e) {
  if (mesh->dim() == 3) return measure_elements_real_tmpl<3>(mesh, a2e);
  if (mesh->dim() == 2) return measure_elements_real_tmpl<2>(mesh, a2e);
  if (mesh->dim() == 1) return measure_elements_real_tmpl<1>(mesh, a2e);
  OMEGA_H_NORETURN(Reals());
}

Reals measure_elements_real(Mesh* mesh) {
  return measure_elements_real(mesh, LOs(mesh->nelems(), 0, 1));
}

}  // end namespace Omega_h
