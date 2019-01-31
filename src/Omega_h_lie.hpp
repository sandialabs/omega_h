#ifndef OMEGA_H_LIE_HPP
#define OMEGA_H_LIE_HPP

/* Lie as in Lie groups and Lie algebras,
 * not as in the opposite of truth.
 */

#include <Omega_h_eigen.hpp>
#include <Omega_h_svd.hpp>

namespace Omega_h {

/* Markley, F. Landis.
   "Unit quaternion from rotation matrix."
   Journal of guidance, control, and dynamics 31.2 (2008): 440-442.

   Modified Shepperd's algorithm to handle input
   tensors that may not be exactly orthogonal */
OMEGA_H_INLINE Vector<4> quaternion_from_tensor(
    Matrix<3, 3> const R) OMEGA_H_NOEXCEPT {
  auto const trR = trace(R);
  auto maxm = trR;
  int maxi = 3;
  Vector<4> q;
  for (int i = 0; i < 3; ++i) {
    if (R(i, i) > maxm) {
      maxm = R(i, i);
      maxi = i;
    }
  }
  if (maxi == 0) {
    q[1] = 1.0 + R(0, 0) - R(1, 1) - R(2, 2);
    q[2] = R(0, 1) + R(1, 0);
    q[3] = R(0, 2) + R(2, 0);
    q[0] = R(2, 1) - R(1, 2);
  } else if (maxi == 1) {
    q[1] = R(1, 0) + R(0, 1);
    q[2] = 1.0 + R(1, 1) - R(2, 2) - R(0, 0);
    q[3] = R(1, 2) + R(2, 1);
    q[0] = R(0, 2) - R(2, 0);
  } else if (maxi == 2) {
    q[1] = R(2, 0) + R(0, 2);
    q[2] = R(2, 1) + R(1, 2);
    q[3] = 1.0 + R(2, 2) - R(0, 0) - R(1, 1);
    q[0] = R(1, 0) - R(0, 1);
  } else if (maxi == 3) {
    q[1] = R(2, 1) - R(1, 2);
    q[2] = R(0, 2) - R(2, 0);
    q[3] = R(1, 0) - R(0, 1);
    q[0] = 1.0 + trR;
  }
  q = normalize(q);
  return q;
}

OMEGA_H_INLINE Vector<3> axis_angle_from_quaternion(
    Vector<4> const q) OMEGA_H_NOEXCEPT {
  auto const divisor = std::sqrt(1.0 - square(q(0)));
  if (divisor < DBL_EPSILON) {
    return zero_vector<3>();
  } else {
    auto const factor = 2.0 * std::acos(q(0)) / divisor;
    Vector<3> aa;
    aa(0) = q(1) * factor;
    aa(1) = q(2) * factor;
    aa(2) = q(3) * factor;
    return aa;
  }
}

// logarithm of a rotation tensor in Special Orthogonal Group(3), as the
// the axis of rotation times the angle of rotation.
OMEGA_H_INLINE Vector<3> axis_angle_from_tensor(
    Matrix<3, 3> const R) OMEGA_H_NOEXCEPT {
  return axis_angle_from_quaternion(quaternion_from_tensor(R));
}

OMEGA_H_INLINE Vector<1> axis_angle_from_tensor(
    Matrix<2, 2> const R) OMEGA_H_NOEXCEPT {
  auto const theta = rotation_angle(R);
  return vector_1(theta);
}

// This function maps a rotation pseudo-vector, aa, to a quaternion, qq
// = (qs, qv), where qv is a vector and qs a scalar, defined as
// follows:
//
//   qv = sin(|aa| / 2) * aa / |aa|
//   qs = cos(|aa| / 2)
//
OMEGA_H_INLINE Vector<4> quaternion_from_axis_angle(
    Vector<3> const aa) OMEGA_H_NOEXCEPT {
  auto const halfnorm = 0.5 * norm(aa);
  auto const temp = 0.5 * sin_x_over_x(halfnorm);
  Vector<4> qq;
  auto const qv = temp * aa;
  for (int i = 0; i < 3; ++i) qq[i + 1] = qv[i];
  qq[0] = std::cos(halfnorm);
  return qq;
}

OMEGA_H_INLINE Matrix<3, 3> tensor_from_quaternion(Vector<4> const qq) {
  auto const qs = qq[0];
  auto const qv = vector_3(qq[1], qq[2], qq[3]);
  auto const I = identity_matrix<3, 3>();
  auto const R = 2.0 * outer_product(qv, qv) + 2.0 * qs * cross(qv) +
                 (2.0 * square(qs) - 1.0) * I;
  return R;
}

OMEGA_H_INLINE Matrix<3, 3> tensor_from_axis_angle(
    Vector<3> const aa) OMEGA_H_NOEXCEPT {
  return tensor_from_quaternion(quaternion_from_axis_angle(aa));
}

OMEGA_H_INLINE Matrix<2, 2> tensor_from_axis_angle(
    Vector<1> const aa) OMEGA_H_NOEXCEPT {
  auto const theta = aa[0];
  return rotate(theta);
}

OMEGA_H_INLINE Matrix<3, 3> pack_polar(
    Matrix<3, 3> const spd, Vector<3> const aa) OMEGA_H_NOEXCEPT {
  Matrix<3, 3> packed = spd;
  packed(0, 1) = aa[0];
  packed(0, 2) = aa[1];
  packed(1, 2) = aa[2];
  return packed;
}

OMEGA_H_INLINE Matrix<3, 3> unpack_polar_spd(
    Matrix<3, 3> const packed) OMEGA_H_NOEXCEPT {
  Matrix<3, 3> spd = packed;
  spd(0, 1) = spd(1, 0);
  spd(0, 2) = spd(2, 0);
  spd(1, 2) = spd(2, 1);
  return spd;
}

OMEGA_H_INLINE Vector<3> unpack_polar_axis_angle(
    Matrix<3, 3> const packed) OMEGA_H_NOEXCEPT {
  Vector<3> aa;
  aa[0] = packed(0, 1);
  aa[1] = packed(0, 2);
  aa[2] = packed(1, 2);
  return aa;
}

OMEGA_H_INLINE Matrix<2, 2> pack_polar(
    Matrix<2, 2> const spd, Vector<1> const aa) OMEGA_H_NOEXCEPT {
  Matrix<2, 2> packed = spd;
  packed(0, 1) = aa[0];
  return packed;
}

OMEGA_H_INLINE Matrix<2, 2> unpack_polar_spd(
    Matrix<2, 2> const packed) OMEGA_H_NOEXCEPT {
  Matrix<2, 2> spd = packed;
  spd(0, 1) = spd(1, 0);
  return spd;
}

OMEGA_H_INLINE Vector<1> unpack_polar_axis_angle(
    Matrix<2, 2> const packed) OMEGA_H_NOEXCEPT {
  Vector<1> aa;
  aa[0] = packed(0, 1);
  return aa;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> log_polar(
    Matrix<dim, dim> const F) OMEGA_H_NOEXCEPT {
  auto const svd = decompose_svd(F);
  auto const X = svd.U;
  auto const D = svd.S;
  auto const YT = svd.V;
  auto const Y = transpose(YT);
  auto const R = X * YT;
  auto const r = axis_angle_from_tensor(R);
  Matrix<dim, dim> d = zero_matrix<dim, dim>();
  for (Int i = 0; i < dim; ++i) d(i, i) = std::log(D(i, i));
  auto const u = Y * d * YT;
  return pack_polar(u, r);
}

OMEGA_H_INLINE Matrix<1, 1> log_polar(Matrix<1, 1> const F) OMEGA_H_NOEXCEPT {
  return matrix_1x1(std::log(F(0, 0)));
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> exp_polar(
    Matrix<dim, dim> const packed_polar) OMEGA_H_NOEXCEPT {
  auto const r = unpack_polar_axis_angle(packed_polar);
  auto const R = tensor_from_axis_angle(r);
  auto const u = unpack_polar_spd(packed_polar);
  auto const U = exp_spd(u);
  auto const F = R * U;
  return F;
}

OMEGA_H_INLINE Matrix<1, 1> exp_polar(
    Matrix<1, 1> const log_F) OMEGA_H_NOEXCEPT {
  return matrix_1x1(std::exp(log_F(0, 0)));
}

// this function modifies the axis-angle rotations of a given array of
// logarithms of polar decompositions such that there are no axis-angle vectors
// which are nearly opposite with angles close to pi. this later allows weighted
// sums of those vectors to give meaningful results and avoids the catastropic
// cancellation case of having (epsilon - pi) and (-epsilon + pi) average to
// zero. tolerance: between 0 and 1, tolerance for opposite pseudo-vectors
// mapping to rotations close to pi
template <Int dim>
OMEGA_H_INLINE void align_packed_axis_angles(Matrix<dim, dim>* const a,
    Int const n, Real const tolerance) OMEGA_H_NOEXCEPT {
  auto const alpha = 1.0 - tolerance;
  auto const s_1 = unpack_polar_axis_angle(a[0]);
  auto s = norm(s_1);
  auto const pi_sq = square(PI);
  auto const two_pi = 2.0 * PI;
  for (Int i = 1; i < n; ++i) {
    auto s_i = unpack_polar_axis_angle(a[i]);
    if ((s_i * s_1) < (-alpha * pi_sq)) {  // pseudo-vectors are nearly opposite
                                           // with angle close to pi
      s_i = s_i - two_pi * normalize(s_i);
      a[i] = pack_polar(unpack_polar_spd(a[i]), s_i);
    }
    s = s + norm(s_i);
  }
  if (s > (n * PI)) {  // renormalize so that ||s_i|| <= pi
    for (Int i = 0; i < n; ++i) {
      auto s_i = unpack_polar_axis_angle(a[i]);
      s_i = s_i - two_pi * normalize(s_i);
      a[i] = pack_polar(unpack_polar_spd(a[i]), s_i);
    }
  }
}

OMEGA_H_INLINE void align_packed_axis_angles(
    Matrix<1, 1>* const, Int const, Real const) OMEGA_H_NOEXCEPT {}

}  // namespace Omega_h

#endif
