#ifndef SIZE_HPP
#define SIZE_HPP

#include "access.hpp"
#include "algebra.hpp"
#include "metric.hpp"
#include "qr.hpp"
#include "space.hpp"

namespace osh {

template <Int sdim, Int edim>
INLINE Few<Vector<sdim>, edim> simplex_basis(Few<Vector<sdim>, edim + 1> p) {
  Few<Vector<sdim>, edim> b;
  for (Int i = 0; i < edim; ++i) b[i] = p[i + 1] - p[0];
  return b;
}

INLINE Real triangle_area(Few<Vector<2>, 2> b) {
  return cross(b[0], b[1]) / 2.0;
}

INLINE Real element_size(Few<Vector<2>, 2> b) { return triangle_area(b); }

INLINE Real triangle_area(Few<Vector<3>, 2> b) {
  return norm(cross(b[0], b[1])) / 2.0;
}

INLINE Real tet_volume(Few<Vector<3>, 3> b) {
  return (cross(b[0], b[1]) * b[2]) / 6.0;
}

INLINE Real element_size(Few<Vector<3>, 3> b) { return tet_volume(b); }

template <Int dim>
INLINE Real iso_edge_length(Few<Vector<dim>, 2> p, Real iso) {
  return norm(p[1] - p[0]) / iso;
}

template <Int dim>
DEVICE Real iso_edge_length(Few<LO, 2> v, Reals coords, Reals isos) {
  auto p = gather_vectors<2, dim>(coords, v);
  auto iso = average(gather_scalars<2>(isos, v));
  return iso_edge_length(p, iso);
}

template <Int dim>
INLINE Real metric_edge_length(Few<Vector<dim>, 2> p, Matrix<dim, dim> metric) {
  return metric_length(metric, p[1] - p[0]);
}

template <Int dim>
DEVICE Real metric_edge_length(Few<LO, 2> v, Reals coords, Reals metrics) {
  auto p = gather_vectors<2, dim>(coords, v);
  auto metric = average_metrics(gather_symms<2, dim>(metrics, v));
  return metric_edge_length(p, metric);
}

template <Int dim>
struct RealEdgeLengths {
  Reals coords;
  RealEdgeLengths(Mesh const* mesh) : coords(mesh->coords()) {}
  DEVICE Real measure(Few<LO, 2> v) const {
    auto p = gather_vectors<2, dim>(coords, v);
    return norm(p[1] - p[0]);
  }
};

template <Int dim>
struct IsoEdgeLengths {
  Reals coords;
  Reals isos;
  IsoEdgeLengths(Mesh const* mesh)
      : coords(mesh->coords()), isos(mesh->get_array<Real>(VERT, "size")) {}
  DEVICE Real measure(Few<LO, 2> v) const {
    return iso_edge_length<dim>(v, coords, isos);
  }
};

template <Int dim>
struct MetricEdgeLengths {
  Reals coords;
  Reals metrics;
  MetricEdgeLengths(Mesh const* mesh)
      : coords(mesh->coords()),
        metrics(mesh->get_array<Real>(VERT, "metric")) {}
  DEVICE Real measure(Few<LO, 2> v) const {
    return metric_edge_length<dim>(v, coords, metrics);
  }
};

Reals measure_edges_real(Mesh* mesh, LOs a2e);
Reals measure_edges_metric(Mesh* mesh, LOs a2e);
Reals measure_edges_real(Mesh* mesh);
Reals measure_edges_metric(Mesh* mesh);

template <Int dim>
INLINE Real real_element_size(Few<Vector<dim>, dim + 1> p) {
  auto b = simplex_basis<dim, dim>(p);
  return element_size(b);
}

struct RealElementSizes {
  Reals coords;
  RealElementSizes(Mesh const* mesh) : coords(mesh->coords()) {}
  template <Int neev>
  DEVICE Real measure(Few<LO, neev> v) const {
    auto p = gather_vectors<neev, neev - 1>(coords, v);
    return real_element_size<neev - 1>(p);
  }
};

Reals measure_elements_real(Mesh* mesh);

INLINE Few<Vector<2>, 3> element_edge_vectors(Few<Vector<2>, 3> p,
                                              Few<Vector<2>, 2> b) {
  Few<Vector<2>, 3> ev;
  ev[0] = b[0];
  ev[1] = p[2] - p[1];
  ev[2] = -b[1];
  return ev;
}

INLINE Few<Vector<3>, 6> element_edge_vectors(Few<Vector<3>, 4> p,
                                              Few<Vector<3>, 3> b) {
  Few<Vector<3>, 6> ev;
  ev[0] = b[0];
  ev[1] = p[2] - p[1];
  ev[2] = -b[1];
  ev[3] = b[2];
  ev[4] = p[3] - p[1];
  ev[5] = p[3] - p[2];
  return ev;
}

INLINE Matrix<2, 2> element_identity_metric(Few<Vector<2>, 3> p) {
  auto b = simplex_basis<2, 2>(p);
  auto ev = element_edge_vectors(p, b);
  Matrix<3, 3> a;
  Vector<3> rhs;
  for (Int i = 0; i < 3; ++i) {
    /* ax^2 + by^2 + 2cxy = 1 */
    a[0][i] = square(ev[i][0]);
    a[1][i] = square(ev[i][1]);
    a[2][i] = 2 * ev[i][0] * ev[i][1];
    rhs[i] = 1.0;
  }
  auto x = invert(a) * rhs;
  return vector2symm(x);
}

INLINE Matrix<3, 3> element_identity_metric(Few<Vector<3>, 4> p) {
  auto b = simplex_basis<3, 3>(p);
  auto ev = element_edge_vectors(p, b);
  Matrix<6, 6> a;
  Vector<6> rhs;
  for (Int i = 0; i < 6; ++i) {
    /* ax^2 + by^2 + cz^2 + 2dxy + 2eyz + 2fxz = 1 */
    a[0][i] = square(ev[i][0]);
    a[1][i] = square(ev[i][1]);
    a[2][i] = square(ev[i][2]);
    a[3][i] = 2 * ev[i][0] * ev[i][1];
    a[4][i] = 2 * ev[i][1] * ev[i][2];
    a[5][i] = 2 * ev[i][0] * ev[i][2];
    rhs[i] = 1.0;
  }
  Vector<6> x;
  /* least squares should decay to exact solution when A is square */
  auto ok = solve_least_squares_qr(a, rhs, x);
  CHECK(ok);
  return vector2symm(x);
}

Reals perfect_size_volumes(Mesh* mesh, Reals v2h);
Reals perfect_metric_volumes(Mesh* mesh, Reals v2m);
Real volume_scalar_for_nelems(Mesh* mesh, Real volume_sum, Real target_nelems);
Reals scale_size_for_nelems(Mesh* mesh, Reals v2h, Real target_nelems);
Reals scale_metric_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems);

}  // end namespace osh

#endif
