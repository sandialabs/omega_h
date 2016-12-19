#ifndef LIE_HPP
#define LIE_HPP

/* Lie as in Lie groups and Lie algebras,
 * not as in the opposite of truth.
 */

#include "eigen.hpp"

namespace Omega_h {

template <Int dim>
INLINE Matrix<dim, dim> log_symm(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::log(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

template <Int dim>
INLINE Matrix<dim, dim> exp_symm(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::exp(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

template <Int dim>
INLINE Matrix<dim, dim> sqrt_symm(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::sqrt(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

INLINE Matrix<3, 3> log_so(Matrix<3, 3> r) {
  auto a = rotation_angle(r);
  if (fabs(a) < EPSILON) return zero_matrix<3,3>();
  if (fabs(a - PI) < EPSILON) {
    auto decomp = decompose_eigen(r);
    auto best_d = fabs(decomp.l[0] - 1.0);
    auto best_i = 0;
    for (Int i = 1; i < 3; ++i) {
      auto d = fabs(decomp.l[i] - 1.0);
      if (d < best_d) {
        best_d = d;
        best_i = i;
      }
    }
    auto v = decomp.q[best_i];
    return PI * cross(v);
  }
  return (a / (2.0 * sin(a))) * (r - transpose(r));
}

INLINE Matrix<3, 3> exp_so(Matrix<3, 3> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::exp(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

template <Int dim>
INLINE Matrix<dim, dim> log_glp(Matrix<dim, dim> a) {
  auto p = sqrt_symm(transpose(a) * a);
  auto u = a * invert(p);
}

template <Int dim>
INLINE Matrix<dim, dim> exp_glp(Matrix<dim, dim> m) {
}

}

#endif
