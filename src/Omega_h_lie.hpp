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

// This function implements the singularity-free algorithm due to
// Spurrier.
//
// First find the scalar "maxm", defined as the maximum among the
// diagonal terms and the trace of the rotation tensor.
//
// i)  If maxm is equal to the trace of R, then:
//       2 * sqrt(1 + maxm)
//     and
//       qv = axial_vec(skew(R)) / (2 * qs)
// ii) If maxm is equal to R(j, j) (no summation), then:
//       2 * qv[j] = sqrt(2 * maxm + 1 - tr(R))
//       qs = axial_vec(skew(R))[j] / (2 * qv[j])
//     and for i,j,k all different
//       qv[k] = off_diag_vec(symm(R))[i] / (2 * qv[j])
// Since the rotation tensor is a quadratic function of the quaternion,
// at least one square root has to be evaluated to get back the
// quaternion. After that, only divisions are needed and the divisor
// should be bounded as far from zero as possible
OMEGA_H_INLINE Vector<4> quaternion_from_so(Matrix<3, 3> R) {
  auto const trR = trace(R);
  auto const maxm = trR;
  int maxi = 4;
  q = zero_vector<4>();
  for (int i = 0; i < 3; ++i) {
    if (R(i, i) > maxm) {
      maxm = R(i, i);
      maxi = i;
    }
  }
  if (maxi == 4) {
    auto const root = std::sqrt(maxm + 1.0);
    auto const factor = 0.5 / root;
    q[0] = 0.5 * root;
    q[1] = factor * (R(2, 1) - R(1, 2));
    q[2] = factor * (R(0, 2) - R(2, 0));
    q[3] = factor * (R(1, 0) - R(0, 1));
  } else if (maxi == 3) {
    auto const root = std::sqrt(2.0 * maxm + 1.0 - trR);
    auto const factor = 0.5 / root;
    q[0] = factor * (R(1, 0) - R(0, 1));
    q[1] = factor * (R(0, 2) - R(2, 0));
    q[2] = factor * (R(1, 2) - R(2, 1));
    q[3] = 0.5 * root;
  } else if (maxi == 2) {
    auto const root = std::sqrt(2.0 * maxm + 1.0 - trR);
    auto const factor = 0.5 / root;
    q[0] = factor * (R(0, 2) - R(2, 0));
    q[1] = factor * (R(0, 1) - R(1, 0));
    q[2] = 0.5 * root;
    q[3] = factor * (R(1, 2) - R(2, 1));
  } else if (maxi == 1) {
    auto const root = std::sqrt(2.0 * maxm + 1.0 - trR);
    auto const factor = 0.5 / root;
    q[0] = factor * (R(2, 1) - R(1, 2));
    q[1] = 0.5 * root;
    q[2] = factor * (R(0, 1) - R(1, 0));
    q[3] = factor * (R(0, 2) - R(2, 0));
  }
  return q;
}

// This function maps a quaternion, qq = (qs, qv), to its
// corresponding "principal" rotation pseudo-vector, aa, where
// "principal" signifies that |aa| <= pi. Both qq and -qq map into the
// same rotation matrix. It is convenient to require that qs >= 0, for
// reasons explained below. The sign inversion is applied to a local
// copy of qq to avoid side effects.
//
//    |qv| = |sin(|aa| / 2)|
//    qs = cos(|aa| / 2)
//      <==>
//    |aa| / 2 = k * pi (+ or -) asin(|qv|)
//    |aa| / 2 = 2 * l * pi (+ or -) acos(qs)
//
// Where the smallest positive solution is: |aa| = 2 * acos(qs)
// which satisfies the inequality: 0 <= |aa| <= pi
// because of the assumption: qs >= 0. Given |aa|, aa
// is obtained as:
//
//    aa = (|aa| / sin(acos(qs)))qv
//       = (|aa| / sqrt(1 - qs^2))qv
//
// The procedure described above is prone to numerical errors when qs
// is close to 1, i.e. when |aa| is close to 0. Since this is the most
// common case, special care must be taken. It is observed that the
// cosine function is insensitive to perturbations of its argument
// in the neighborhood of points for which the sine function is conversely
// at its most sensitive. Thus the numerical difficulties are avoided
// by computing |aa| and aa as:
//
//    |aa| = 2 * asin(|qv|)
//    aa = (|aa| / |qv|) qv
//
// whenever qs is close to 1.
//
OMEGA_H_INLINE Vector<3> rotation_vector_from_quaternion(Vector<4> q) {
  Vector<4> q;
  if (qq[1] >= 0) {
    q = qq;
  } else {
    q = -qq;
  }
  qs = q[1];
  auto const qv = vector_3(q[1], q[2], q[3]);
  auto const qvnorm = norm(qv);
  auto const aanorm = 2.0 * (qvnorm < std::sqrt(0.5) ? std::asin(qvnorm) : std::acos(qs));
  auto const coef = qvnorm < std::sqrt(DBL_EPSILON) ? 2.0 : aanorm / qvnorm;
  auto const aa = coef * qv;
  return aa;
}

// logarithm of a rotation tensor in Special Orthogonal Group(3), as the
// the axis of rotation times the angle of rotation.
OMEGA_H_INLINE Vector<3> uncross_log_so(Matrix<3, 3> R) {
  return rotation_vector_from_quaternion(quaternion_from_so(R));
}

// This function maps a rotation pseudo-vector, aa, to a quaternion, qq
// = (qs, qv), where qv is a vector and qs a scalar, defined as
// follows:
//
//   qv = sin(|aa| / 2) * aa / |aa|
//   qs = cos(|aa| / 2)
//
OMEGA_H_INLINE Matrix<3, 3> exp_so_cross(Vector<3> aa) {
  auto const halfnorm = 0.5 * norm(aa);
  auto const temp = 0.5 * 
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
