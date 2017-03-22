#ifndef OMEGA_H_SIZE_HPP
#define OMEGA_H_SIZE_HPP

#include "Omega_h_metric.hpp"
#include "Omega_h_qr.hpp"
#include "Omega_h_simplex.hpp"
#include "Omega_h_adj.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

template <Int sdim, Int edim>
OMEGA_H_INLINE Matrix<sdim, edim> simplex_basis(Few<Vector<sdim>, edim + 1> p) {
  Matrix<sdim, edim> b;
  for (Int i = 0; i < edim; ++i) b[i] = p[i + 1] - p[0];
  return b;
}

template <Int dim>
struct Affine {
  Matrix<dim, dim> r;
  Vector<dim> t;
};

template <Int dim>
OMEGA_H_INLINE Vector<dim> operator*(Affine<dim> a, Vector<dim> v) {
  return (a.r * v) + a.t;
}

template <Int dim>
OMEGA_H_INLINE Affine<dim> invert(Affine<dim> a) {
  Affine<dim> ai;
  ai.r = invert(a.r);
  ai.t = -(ai.r * a.t);
  return ai;
}

template <Int dim>
OMEGA_H_INLINE Affine<dim> simplex_affine(Few<Vector<dim>, dim + 1> p) {
  Affine<dim> a;
  a.r = simplex_basis<dim, dim>(p);
  a.t = p[0];
  return a;
}

template <Int dim>
struct Sphere {
  Vector<dim> c;
  Real r;
};

template <Int dim>
struct Plane {
  Vector<dim> n; /*!< Unit-length normal vector. */
  Real d;        /*!< Signed perpendicular distance to the origin. */
};

template <Int dim>
OMEGA_H_INLINE Real distance(Plane<dim> plane, Vector<dim> point) {
  return (plane.n * point) - plane.d;
}

template <Int dim>
OMEGA_H_INLINE Vector<dim + 1> form_barycentric(Vector<dim> c) {
  Vector<dim + 1> xi;
  xi[dim] = 1.0;
  for (Int i = 0; i < dim; ++i) {
    xi[i] = c[i];
    xi[dim] -= c[i];
  }
  return xi;
}

template <Int n>
OMEGA_H_INLINE bool is_barycentric_inside(Vector<n> xi) {
  return 0.0 <= minimum(xi) && maximum(xi) <= 1.0;
}

OMEGA_H_INLINE Real triangle_area(Few<Vector<2>, 2> b) {
  return cross(b[0], b[1]) / 2.0;
}

OMEGA_H_INLINE Real element_size(Few<Vector<2>, 2> b) { return triangle_area(b); }

OMEGA_H_INLINE Real triangle_area(Few<Vector<3>, 2> b) {
  return norm(cross(b[0], b[1])) / 2.0;
}

OMEGA_H_INLINE Real tet_volume(Few<Vector<3>, 3> b) {
  return (cross(b[0], b[1]) * b[2]) / 6.0;
}

OMEGA_H_INLINE Real element_size(Few<Vector<3>, 3> b) { return tet_volume(b); }

/* Loseille, Adrien, and Rainald Lohner.
 * "On 3D anisotropic local remeshing for surface, volume and boundary layers."
 * Proceedings of the 18th International Meshing Roundtable.
 * Springer Berlin Heidelberg, 2009. 611-630.
 *
 * Loseille's edge length integral assumes an interpolation $h(t) =
 * h_0^{1-t}h_1^t$,
 * which is consistent with the Log-Euclidean metric interpolation we now use.
 */

OMEGA_H_INLINE Real edge_length(Real l_a, Real l_b) {
  if (::fabs(l_a - l_b) > 1e-3) {
    return (l_a - l_b) / (::log(l_a / l_b));
  }
  return (l_a + l_b) / 2.;
}

template <Int space_dim, Int metric_dim>
OMEGA_H_INLINE Real squared_metric_length(
    Vector<space_dim> v, Matrix<metric_dim, metric_dim> m) {
  return metric_product(m, v);
}

template <Int space_dim>
OMEGA_H_INLINE Real squared_metric_length(Vector<space_dim> v, NoMetric) {
  return norm_squared(v);
}

template <Int space_dim, Int metric_dim>
OMEGA_H_INLINE Real metric_edge_length(
    Few<Vector<space_dim>, 2> p, Few<Matrix<metric_dim, metric_dim>, 2> ms) {
  auto v = p[1] - p[0];
  auto l_a = metric_length(ms[0], v);
  auto l_b = metric_length(ms[1], v);
  return edge_length(l_a, l_b);
}

template <Int space_dim, Int metric_dim>
OMEGA_H_INLINE Real metric_edge_length(Few<LO, 2> v, Reals coords, Reals metrics) {
  auto p = gather_vectors<2, space_dim>(coords, v);
  auto ms = gather_symms<2, metric_dim>(metrics, v);
  return metric_edge_length<space_dim, metric_dim>(p, ms);
}

template <Int space_dim>
struct RealEdgeLengths {
  Reals coords;
  RealEdgeLengths(Mesh const* mesh) : coords(mesh->coords()) {}
  OMEGA_H_DEVICE Real measure(Few<LO, 2> v) const {
    auto p = gather_vectors<2, space_dim>(coords, v);
    return norm(p[1] - p[0]);
  }
};

template <Int space_dim, Int metric_dim>
struct MetricEdgeLengths {
  Reals coords;
  Reals metrics;
  MetricEdgeLengths(Mesh const* mesh)
      : coords(mesh->coords()),
        metrics(mesh->get_array<Real>(VERT, "metric")) {}
  OMEGA_H_DEVICE Real measure(Few<LO, 2> v) const {
    return metric_edge_length<space_dim, metric_dim>(v, coords, metrics);
  }
};

Reals measure_edges_real(Mesh* mesh, LOs a2e);
Reals measure_edges_metric(Mesh* mesh, LOs a2e);
Reals measure_edges_real(Mesh* mesh);
Reals measure_edges_metric(Mesh* mesh);

template <Int dim>
OMEGA_H_INLINE Real real_element_size(Few<Vector<dim>, dim + 1> p) {
  auto b = simplex_basis<dim, dim>(p);
  return element_size(b);
}

struct RealElementSizes {
  Reals coords;
  RealElementSizes(Mesh const* mesh) : coords(mesh->coords()) {}
  template <Int neev>
  OMEGA_H_DEVICE Real measure(Few<LO, neev> v) const {
    auto p = gather_vectors<neev, neev - 1>(coords, v);
    return real_element_size<neev - 1>(p);
  }
};

Reals measure_elements_real(Mesh* mesh);

OMEGA_H_INLINE Few<Vector<2>, 3> element_edge_vectors(
    Few<Vector<2>, 3> p, Few<Vector<2>, 2> b) {
  Few<Vector<2>, 3> ev;
  ev[0] = b[0];
  ev[1] = p[2] - p[1];
  ev[2] = -b[1];
  return ev;
}

OMEGA_H_INLINE Few<Vector<3>, 6> element_edge_vectors(
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

template <typename EdgeVectors, typename Metric>
OMEGA_H_INLINE Real mean_squared_metric_length(
    EdgeVectors edge_vectors, Metric metric) {
  auto nedges = EdgeVectors::size;
  Real msl = 0;
  for (Int i = 0; i < nedges; ++i) {
    msl += squared_metric_length(edge_vectors[i], metric);
  }
  return msl / nedges;
}

template <typename EdgeVectors>
OMEGA_H_INLINE Real mean_squared_real_length(EdgeVectors edge_vectors) {
  return mean_squared_metric_length(edge_vectors, NoMetric());
}

template <Int dim>
OMEGA_H_INLINE Real element_implied_length(Few<Vector<dim>, dim + 1> p) {
  auto b = simplex_basis<dim, dim>(p);
  auto ev = element_edge_vectors(p, b);
  auto h = sqrt(mean_squared_real_length(ev));
  return h;
}

OMEGA_H_INLINE Matrix<2, 2> element_implied_metric(Few<Vector<2>, 3> p) {
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

OMEGA_H_INLINE Matrix<3, 3> element_implied_metric(Few<Vector<3>, 4> p) {
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

OMEGA_H_INLINE Vector<1> get_side_normal(Few<Vector<1>, 2>, Int ivert) {
  return vector_1((ivert == 1) ? 1.0 : -1.0);
}

OMEGA_H_INLINE Vector<2> get_side_normal(Few<Vector<2>, 3> p, Int iedge) {
  auto a = p[down_template(2, 1, iedge, 0)];
  auto b = p[down_template(2, 1, iedge, 1)];
  return -perp(b - a);
}

OMEGA_H_INLINE Vector<3> get_side_normal(Few<Vector<3>, 4> p, Int iface) {
  auto a = p[down_template(3, 2, iface, 0)];
  auto b = p[down_template(3, 2, iface, 1)];
  auto c = p[down_template(3, 2, iface, 2)];
  return cross(b - a, c - a);
}

template <Int dim>
OMEGA_H_INLINE Plane<dim> get_side_plane(Few<Vector<dim>, dim + 1> p, Int iside) {
  auto n = get_side_normal(p, iside);
  auto a = p[down_template(dim, dim - 1, iside, 0)];
  return {n, n * a};
}

template <Int dim>
OMEGA_H_INLINE Sphere<dim> get_insphere(Few<Vector<dim>, dim + 1> p) {
  auto nsides = dim + 1;
  Few<Plane<dim>, dim + 1> planes;
  for (Int iside = 0; iside < nsides; ++iside) {
    planes[iside] = get_side_plane(p, iside);
  }
  Matrix<dim, dim> a;
  Vector<dim> b;
  for (Int i = 0; i < dim; ++i) {
    a[i] = planes[dim].n - planes[i].n;
    b[i] = planes[dim].d - planes[i].d;
  }
  auto c = invert(transpose(a)) * b;
  auto r = -distance(planes[0], c) / norm(planes[0].n);
  return {c, r};
}

template <>
OMEGA_H_INLINE Sphere<1> get_insphere(Few<Vector<1>, 2> p) {
  return {average(p), fabs((p[1] - p[0])[0] / 2.0)};
}

template <Int dim>
Vector<dim> get_volume_vert_gradient(Few<Vector<dim>, dim + 1> p, Int ivert) {
  auto iside = opposite_template(dim, VERT, ivert);
  auto n = -get_side_normal(p, iside);
  return n / Real(factorial(dim));
}

/* This code is copied from the tricircumcenter3d() function
 * by Shewchuk:
 * http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
 * To:             compgeom-discuss@research.bell-labs.com
 * Subject:        Re: circumsphere
 * Date:           Wed, 1 Apr 98 0:34:28 EST
 * From:           Jonathan R Shewchuk <jrs+@cs.cmu.edu>
 *
 * given the basis vectors of a triangle in 3D,
 * this function returns the vector from the first vertex
 * to the triangle's circumcenter
 */

OMEGA_H_INLINE Vector<3> get_circumcenter_vector(Few<Vector<3>, 2> basis) {
  auto ba = basis[0];
  auto ca = basis[1];
  auto balength = norm_squared(ba);
  auto calength = norm_squared(ca);
  auto crossbc = cross(ba, ca);
  auto factor = 0.5 / norm_squared(crossbc);
  auto xcirca = ((balength * ca[1] - calength * ba[1]) * crossbc[2] -
                    (balength * ca[2] - calength * ba[2]) * crossbc[1]) *
                factor;
  auto ycirca = ((balength * ca[2] - calength * ba[2]) * crossbc[0] -
                    (balength * ca[0] - calength * ba[0]) * crossbc[2]) *
                factor;
  auto zcirca = ((balength * ca[0] - calength * ba[0]) * crossbc[1] -
                    (balength * ca[1] - calength * ba[1]) * crossbc[0]) *
                factor;
  return vector_3(xcirca, ycirca, zcirca);
}

template <Int dim>
OMEGA_H_INLINE Vector<dim> get_triangle_normal(
    Vector<dim> a, Vector<dim> b, Vector<dim> c) {
  return cross(b - a, c - a);
}

}  // end namespace Omega_h

#endif
