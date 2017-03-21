#ifndef LIE_HPP
#define LIE_HPP

/* Lie as in Lie groups and Lie algebras,
 * not as in the opposite of truth.
 */

#include "eigen.hpp"

namespace Omega_h {

// logarithm of a symmetric positive definite tensor
template <Int dim>
INLINE Matrix<dim, dim> log_spd(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::log(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// exponential resulting in a symmetric positive definite tensor
template <Int dim>
INLINE Matrix<dim, dim> exp_spd(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::exp(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// square root of a symmetric positive definite tensor
template <Int dim>
INLINE Matrix<dim, dim> sqrt_spd(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = ::sqrt(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// logarithm of a tensor in Special Orthogonal Group(3), as axis times angle
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
  return (a / sin(a)) * uncross(r - transpose(r));
}

// exponential of axis-angle, resulting in an SO(3) tensor
INLINE Matrix<3, 3> exp_so(Vector<3> axis_angle) {
  auto a = norm(axis_angle);
  if (fabs(a) < EPSILON) return identity_matrix<3, 3>();
  return rotate(a, axis_angle / a);
}

// logarithm of a tensor in Special Orthogonal Group(2), as angle
INLINE Real log_so(Matrix<2, 2> r) { return rotation_angle(r); }

// exponential of angle, resulting in an SO(2) tensor
INLINE Matrix<2, 2> exp_so(Real angle) { return rotate(angle); }

template <Int dim>
struct LogDecomp;

template <>
struct LogDecomp<2> {
  Real log_u;
  Matrix<2, 2> log_p;
};

template <>
struct LogDecomp<3> {
  Vector<3> log_u;
  Matrix<3, 3> log_p;
};

OMEGA_H_INLINE constexpr Int polar_dofs(Int dim) { return matrix_dofs(dim); }

template <Int dim>
struct LogGLP;

template <Int dim>
INLINE LogDecomp<dim> log_glp(Matrix<dim, dim> a) {
  auto aa_dc = decompose_eigen(transpose(a) * a);
  Vector<dim> p_l;
  for (Int i = 0; i < dim; ++i) p_l[i] = ::sqrt(aa_dc.l[i]);
  auto p = compose_eigen(aa_dc.q, p_l);
  auto u = a * invert(p);
  auto log_u = log_so(u);
  Vector<dim> log_p_l;
  for (Int i = 0; i < dim; ++i) log_p_l[i] = ::log(p_l[i]);
  auto log_p = compose_ortho(aa_dc.q, log_p_l);
  return {log_u, log_p};
}

template <Int dim>
INLINE Matrix<dim, dim> exp_glp(LogDecomp<dim> log_a) {
  auto p = exp_spd(log_a.log_p);
  auto u = exp_so(log_a.log_u);
  return u * p;
}
}  // namespace Omega_h

#endif
