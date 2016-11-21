#ifndef SIZE_HPP
#define SIZE_HPP

#include "access.hpp"
#include "algebra.hpp"
#include "metric.hpp"
#include "qr.hpp"
#include "space.hpp"

namespace Omega_h {

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

/* Loseille, Adrien, and Rainald Lohner.
 * "On 3D anisotropic local remeshing for surface, volume and boundary layers."
 * Proceedings of the 18th International Meshing Roundtable.
 * Springer Berlin Heidelberg, 2009. 611-630.
 *
 * Loseille's edge length integral assumes an interpolation $h(t) =
 * h_0^{1-t}h_1^t$,
 * which is consistent with the Log-Euclidean metric interpolation we now use.
 */

INLINE Real edge_length(Real l_a, Real l_b) {
  if (::fabs(l_a - l_b) > 1e-3) {
    return (l_a - l_b) / (::log(l_a / l_b));
  }
  return (l_a + l_b) / 2.;
}

template <Int dim>
INLINE Real iso_edge_length(Few<Vector<dim>, 2> p, Few<Real, 2> hs) {
  auto real_l = norm(p[1] - p[0]);
  auto l_a = real_l / hs[0];
  auto l_b = real_l / hs[0];
  return edge_length(l_a, l_b);
}

template <Int dim>
DEVICE Real iso_edge_length(Few<LO, 2> v, Reals coords, Reals isos) {
  auto p = gather_vectors<2, dim>(coords, v);
  auto hs = gather_scalars<2>(isos, v);
  return iso_edge_length(p, hs);
}

template <Int dim>
INLINE Real metric_edge_length(
    Few<Vector<dim>, 2> p, Few<Matrix<dim, dim>, 2> ms) {
  auto v = p[1] - p[0];
  auto l_a = metric_length(ms[0], v);
  auto l_b = metric_length(ms[1], v);
  return edge_length(l_a, l_b);
}

template <Int dim>
DEVICE Real metric_edge_length(Few<LO, 2> v, Reals coords, Reals metrics) {
  auto p = gather_vectors<2, dim>(coords, v);
  auto ms = gather_symms<2, dim>(metrics, v);
  return metric_edge_length(p, ms);
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

INLINE Few<Vector<2>, 3> element_edge_vectors(
    Few<Vector<2>, 3> p, Few<Vector<2>, 2> b) {
  Few<Vector<2>, 3> ev;
  ev[0] = b[0];
  ev[1] = p[2] - p[1];
  ev[2] = -b[1];
  return ev;
}

INLINE Few<Vector<3>, 6> element_edge_vectors(
    Few<Vector<3>, 4> p, Few<Vector<3>, 3> b) {
  Few<Vector<3>, 6> ev;
  ev[0] = b[0];
  ev[1] = p[2] - p[1];
  ev[2] = -b[1];
  ev[3] = b[2];
  ev[4] = p[3] - p[1];
  ev[5] = p[3] - p[2];
  return ev;
}

template <Int dim>
INLINE Real squared_metric_length(Vector<dim> v, DummyIsoMetric) {
  return norm_squared(v);
}

template <Int dim>
INLINE Real squared_metric_length(Vector<dim> v, Matrix<dim, dim> m) {
  return metric_product(m, v);
}

template <typename EdgeVectors, typename Metric>
INLINE Real mean_squared_metric_length(
    EdgeVectors edge_vectors, Metric metric) {
  auto nedges = EdgeVectors::size;
  Real msl = 0;
  for (Int i = 0; i < nedges; ++i) {
    msl += squared_metric_length(edge_vectors[i], metric);
  }
  return msl / nedges;
}

template <typename EdgeVectors>
INLINE Real mean_squared_real_length(EdgeVectors edge_vectors) {
  return mean_squared_metric_length(edge_vectors, DummyIsoMetric());
}

template <Int dim>
INLINE Real element_implied_size(Few<Vector<dim>, dim + 1> p) {
  auto b = simplex_basis<dim, dim>(p);
  auto ev = element_edge_vectors(p, b);
  auto h = sqrt(mean_squared_real_length(ev));
  return h;
}

INLINE Matrix<2, 2> element_implied_metric(Few<Vector<2>, 3> p) {
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

INLINE Matrix<3, 3> element_implied_metric(Few<Vector<3>, 4> p) {
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
  auto x = solve_using_qr(a, rhs);
  return vector2symm(x);
}

template <Int dim>
struct ParentElementSize;

template <>
struct ParentElementSize<2> {
  static constexpr Real value = 1.0 / 2.0;
};

template <>
struct ParentElementSize<3> {
  static constexpr Real value = 1.0 / 6.0;
};

Reals get_mident_isos(Mesh* mesh, Int ent_dim, LOs entities, Reals v2h);
Reals interpolate_between_isos(Reals a, Reals b, Real t);
Reals linearize_isos(Reals isos);
Reals delinearize_isos(Reals log_isos);
Reals project_isos(Mesh* mesh, Reals e2h);

}  // end namespace Omega_h

#endif
