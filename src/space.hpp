#ifndef SPACE_HPP
#define SPACE_HPP

#include "algebra.hpp"

namespace Omega_h {

INLINE Real cross(Vector<2> a, Vector<2> b) {
  return (a[0] * b[1] - a[1] * b[0]);
}

INLINE Vector<3> uncross(Matrix<3, 3> c) {
  return vector_3(c[1][2] - c[2][1], c[2][0] - c[0][2], c[0][1] - c[1][0]) / 2.;
}

/* Rodrigues' Rotation Formula */
INLINE Matrix<3, 3> rotate(Real angle, Vector<3> axis) {
  return cos(angle) * identity_matrix<3, 3>() + sin(angle) * cross(axis) +
         (1 - cos(angle)) * tensor_product(axis, axis);
}

INLINE Matrix<2, 2> rotate(Real angle) {
  return matrix_2x2(cos(angle), -sin(angle), sin(angle), cos(angle));
}

INLINE Real rotation_angle(Matrix<2, 2> r) { return acos(r[0][0]); }

INLINE Real rotation_angle(Matrix<3, 3> r) __attribute__((pure));
INLINE Real rotation_angle(Matrix<3, 3> r) {
  return acos((trace(r) - 1.0) / 2.0);
}

INLINE Vector<2> perp(Vector<2> v) { return vector_2(-v[1], v[0]); }

INLINE Real determinant(Matrix<2, 2> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return a * d - b * c;
}

INLINE Real determinant(Matrix<3, 3> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[2][0];
  Real d = m[0][1];
  Real e = m[1][1];
  Real f = m[2][1];
  Real g = m[0][2];
  Real h = m[1][2];
  Real i = m[2][2];
  return (a * e * i) + (b * f * g) + (c * d * h) - (c * e * g) - (b * d * i) -
         (a * f * h);
}

INLINE Matrix<2, 2> invert(Matrix<2, 2> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return matrix_2x2(d, -b, -c, a) / determinant(m);
}

INLINE Matrix<3, 3> invert(Matrix<3, 3> a) {
  Matrix<3, 3> b;
  b[0] = cross(a[1], a[2]);
  b[1] = cross(a[2], a[0]);
  b[2] = cross(a[0], a[1]);
  return transpose(b) / determinant(a);
}

INLINE Matrix<2, 2> form_ortho_basis(Vector<2> v) {
  Matrix<2, 2> A;
  A[0] = normalize(v);
  A[1] = perp(A[0]);
  return A;
}

INLINE Matrix<3, 3> form_ortho_basis(Vector<3> v) {
  Matrix<3, 3> A;
  A[0] = v;
  /* tiny custom code to sort components by absolute value */
  struct {
    Int i;
    Real m;
  } s[3] = {{0, fabs(v[0])}, {1, fabs(v[1])}, {2, fabs(v[2])}};
  if (s[2].m > s[1].m) swap2(s[1], s[2]);
  if (s[1].m > s[0].m) swap2(s[0], s[1]);
  if (s[2].m > s[1].m) swap2(s[1], s[2]);
  /* done, components sorted by increasing magnitude */
  Int lc = s[0].i;
  Int mc = s[1].i;
  Int sc = s[2].i;
  /* use the 2D rotation on the largest components
     (rotate v around the smallest axis) */
  A[1][lc] = -v[mc];
  A[1][mc] = v[lc];
  /* and make the last component zero so that A[0] * A[1] == 0 */
  A[1][sc] = 0;
  A[1] = normalize(A[1]);
  /* now we have 2 orthogonal unit vectors, cross product gives the third */
  A[2] = cross(A[0], A[1]);
  return A;
}

template <Int dim>
struct Affine {
  Matrix<dim, dim> r;
  Vector<dim> t;
};

template <Int dim>
INLINE Vector<dim> operator*(Affine<dim> a, Vector<dim> v) {
  return (a.r * v) + a.t;
}

template <Int dim>
INLINE Affine<dim> invert(Affine<dim> a) {
  Affine<dim> ai;
  ai.r = invert(a.r);
  ai.t = -(ai.r * a.t);
  return ai;
}

template <Int sdim, Int edim>
INLINE Matrix<sdim, edim> simplex_basis(Few<Vector<sdim>, edim + 1> p) {
  Matrix<sdim, edim> b;
  for (Int i = 0; i < edim; ++i) b[i] = p[i + 1] - p[0];
  return b;
}

template <Int dim>
INLINE Affine<dim> simplex_affine(Few<Vector<dim>, dim + 1> p) {
  Affine<dim> a;
  a.r = simplex_basis<dim, dim>(p);
  a.t = p[0];
  return a;
}

template <Int dim>
INLINE Vector<dim + 1> form_barycentric(Vector<dim> c) {
  Vector<dim + 1> xi;
  xi[dim] = 1.0;
  for (Int i = 0; i < dim; ++i) {
    xi[i] = c[i];
    xi[dim] -= c[i];
  }
  return xi;
}

template <Int n>
INLINE bool is_barycentric_inside(Vector<n> xi) {
  return 0.0 <= minimum(xi) && maximum(xi) <= 1.0;
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

INLINE Vector<3> get_circumcenter_vector(Few<Vector<3>, 2> basis) {
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

}  // end namespace Omega_h

#endif
