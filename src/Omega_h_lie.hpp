#ifndef OMEGA_H_LIE_HPP
#define OMEGA_H_LIE_HPP

/* Lie as in Lie groups and Lie algebras,
 * not as in the opposite of truth.
 */

#include <Omega_h_eigen.hpp>

namespace Omega_h {

// logarithm of a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> log_spd(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen_jacobi(m);
  for (Int i = 0; i < dim; ++i) decomp.l[i] = std::log(decomp.l[i]);
  return compose_ortho(decomp.q, decomp.l);
}

// exponential resulting in a symmetric positive definite tensor
template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> exp_spd(Matrix<dim, dim> m) {
  auto decomp = decompose_eigen_jacobi(m);
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
  // R = cos(theta) * I + sin(theta) * cross(v) + (1 - cos(theta)) * outer_product(u, u)
  // R - R^T = 2 * sin(theta) * cross(v)
  return (theta / std::sin(theta)) * (r - transpose(r));
}

// exponential of axis-angle, resulting in an SO(3) tensor
OMEGA_H_INLINE Matrix<3, 3> exp_so(Matrix<3> log_r) {
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

// exponential of angle, resulting in an SO(2) tensor
OMEGA_H_INLINE Matrix<2, 2> exp_so(Matrix<2, 2> log_r) {
  auto const theta = 0.5 * (log_r(1, 0) - log_r(0, 1));
  return rotate(theta);
}

/* get the logarithm of a tensor in the "identity component of the general linear group",
   denoted by GL+(n) in: 

   Mota, Alejandro, et al.
   "Lie-group interpolation and variational recovery for internal variables."
   Computational Mechanics 52.6 (2013): 1281-1299.

   The tensor A is first polar-decomposed into a rotation R and a symmetric tensor S
   such that A=RS.
   Then the logarithms of the rotation and symmetric tensor are taken separately.

   Finally, we note that the logarithm of a rotation in SO(n) has ((n - 1) * n / 2)
 */

template <Int dim>
OMEGA_H_INLINE Vector<dim * dim> log_glp(Matrix<dim, dim> A) {
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
}

}  // namespace Omega_h

#endif
