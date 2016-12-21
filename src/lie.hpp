#ifndef LIE_HPP
#define LIE_HPP

/* Lie as in Lie groups and Lie algebras,
 * not as in the opposite of truth.
 */

#include "eigen.hpp"
#include "access.hpp"

namespace Omega_h {

// logarithm of a symmetrix positive definite tensor
template <Int dim>
INLINE Matrix<dim, dim> log_spd(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::log(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// exponential of a symmetrix positive definite tensor
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

INLINE Vector<3> log_so(Matrix<3, 3> r) {
  auto a = rotation_angle(r);
  if (fabs(a) < EPSILON) return zero_vector<3>();
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
    return PI * v;
  }
  return a * (uncross(r - transpose(r)) / sin(a));
}

INLINE Matrix<3, 3> exp_so(Matrix<3, 3> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::exp(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

template <Int dim>
struct PolarDecomp {
  Matrix<dim,dim> u; // the unitary (rotational) tensor
  Matrix<dim,dim> p; // the Symmetric Positive Definite tensor
};

OMEGA_H_INLINE constexpr Int polar_dofs(Int dim) {
  return matrix_dofs(sim) + symm_dofs(dim);
}

template <Int dim>
INLINE PolarDecomp<dim> log_glp(Matrix<dim, dim> a) {
  auto p = sqrt_symm(transpose(a) * a);
  auto u = a * invert(p);
  return {u, p};
}

template <Int dim>
INLINE Matrix<dim, dim> exp_glp(PolarDecomp<dim> log_a) {
}

}

#endif
