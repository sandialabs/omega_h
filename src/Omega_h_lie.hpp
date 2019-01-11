#ifndef OMEGA_H_LIE_HPP
#define OMEGA_H_LIE_HPP

/* Lie as in Lie groups and Lie algebras,
 * not as in the opposite of truth.
 */

#include <Omega_h_eigen.hpp>
#include <Omega_h_svd.hpp>

namespace Omega_h {

// logarithm of a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE_BIG Matrix<dim, dim> log_spd_old(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::log(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// exponential resulting in a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE_BIG Matrix<dim, dim> exp_spd_old(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::exp(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// square root of a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> sqrt_spd(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::sqrt(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

/* logarithm of a tensor in Special Orthogonal Group(3), as the
   skew-symmetric cross product tensor of the axis of rotation times
   the angle of rotation.
  */
OMEGA_H_INLINE Matrix<3, 3> log_so(Matrix<3, 3> r) {
  auto theta = rotation_angle(r);
  if (std::abs(theta) < EPSILON) return zero_matrix<3, 3>();
  if (std::abs(theta - PI) < EPSILON) {
    auto decomp = decompose_eigen_jacobi(r);
    auto best_d = std::abs(decomp.l[0] - 1.0);
    auto best_i = 0;
    for (Int i = 1; i < 3; ++i) {
      auto d = std::abs(decomp.l[i] - 1.0);
      if (d < best_d) {
        best_d = d;
        best_i = i;
      }
    }
    auto v = decomp.q[best_i];
    return PI * cross(v);
  }
  // R = cos(theta) * I + sin(theta) * cross(v) + (1 - cos(theta)) *
  // outer_product(u, u)
  // R - R^T = 2 * sin(theta) * cross(v)
  return (theta / (2.0 * std::sin(theta))) * (r - transpose(r));
}

// exponential of axis-angle, resulting in an SO(3) tensor
OMEGA_H_INLINE Matrix<3, 3> exp_so(Matrix<3, 3> log_r) {
  auto const v_times_theta = uncross(log_r);
  auto const theta = norm(v_times_theta);
  if (std::abs(theta) < EPSILON) return identity_matrix<3, 3>();
  auto const v = v_times_theta * (1.0 / theta);
  return rotate(theta, v);
}

// logarithm of a tensor in Special Orthogonal Group(2)
OMEGA_H_INLINE Matrix<2, 2> log_so(Matrix<2, 2> r) {
  auto const theta = rotation_angle(r);
  return matrix_2x2(0, -theta, theta, 0);
}

// exponential resulting in an SO(2) tensor
OMEGA_H_INLINE Matrix<2, 2> exp_so(Matrix<2, 2> log_r) {
  auto const theta = 0.5 * (log_r(1, 0) - log_r(0, 1));
  return rotate(theta);
}

OMEGA_H_INLINE Matrix<1, 1> log_so(Matrix<1, 1>) { return matrix_1x1(0.0); }

OMEGA_H_INLINE Matrix<1, 1> exp_so(Matrix<1, 1>) { return matrix_1x1(1.0); }

/* get the logarithm of a tensor in the "identity component of the general
   linear group",
   denoted by GL+(n) in:

   Mota, Alejandro, et al.
   "Lie-group interpolation and variational recovery for internal variables."
   Computational Mechanics 52.6 (2013): 1281-1299.

   We deviate from the approach outlined in the above paper.
   Instead of dealing with R and S, we take logarithms of the three SVD
   components
   U, D, and V.
 */

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> log_glp(Matrix<dim, dim> A) {
  // A = U * D * V^T
  // A * A^T = U * D^2 * U^T
  auto const UD = decompose_eigen_jacobi(A * transpose(A));
  auto const U = UD.q;
  auto const D_sq = UD.l;
  Vector<dim> D, D_inv, log_D;
  for (Int i = 0; i < dim; ++i) {
    D[i] = std::sqrt(D_sq[i]);
    D_inv[i] = 1.0 / D[i];
    log_D[i] = std::log(D[i]);
  }
  // V^T = D^{-1} * U^T * A
  auto const VT = diagonal(D_inv) * transpose(U) * A;
  auto const V = transpose(VT);
  auto const log_U = log_so(U);
  auto const log_V = log_so(V);
  Matrix<dim, dim> log_A;
  // mix the independent components as follows:
  // the lower triangle will store the lower triangle
  // of log(U), the upper triangle will be upper triangle of log(V),
  // and the diagonal will be the diagonal of log(D)
  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      if (i < j)
        log_A(i, j) = log_V(i, j);
      else if (i > j)
        log_A(i, j) = log_U(i, j);
      else
        log_A(i, i) = log_D(i);
    }
  }
  return log_A;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> exp_glp(Matrix<dim, dim> log_A) {
  auto log_U = zero_matrix<dim, dim>();
  auto log_V = zero_matrix<dim, dim>();
  Vector<dim> log_D;
  for (Int i = 0; i < dim; ++i) {
    for (Int j = 0; j < dim; ++j) {
      if (i < j)
        log_V(i, j) = log_A(i, j);
      else if (i > j)
        log_U(i, j) = log_A(i, j);
      else
        log_D(i) = log_A(i, i);
    }
  }
  log_U -= transpose(log_U);
  log_V -= transpose(log_V);
  auto const U = exp_so(log_U);
  auto const V = exp_so(log_V);
  Vector<dim> D;
  for (Int i = 0; i < dim; ++i) D(i) = std::exp(log_D(i));
  return U * diagonal(D) * transpose(V);
}

}  // namespace Omega_h

#endif
