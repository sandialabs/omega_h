#ifndef OMEGA_H_LIE_HPP
#define OMEGA_H_LIE_HPP

/* Lie as in Lie groups and Lie algebras,
 * not as in the opposite of truth.
 */

#include <Omega_h_eigen.hpp>
#include <Omega_h_svd.hpp>

namespace Omega_h {

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
OMEGA_H_INLINE Vector<4> quaternion_from_tensor(Matrix<3, 3> const R) OMEGA_H_NOEXCEPT {
  auto const trR = trace(R);
  auto maxm = trR;
  int maxi = 4;
  Vector<4> q;
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
OMEGA_H_INLINE Vector<3> axis_angle_from_quaternion(Vector<4> const qq) OMEGA_H_NOEXCEPT {
  Vector<4> q;
  if (qq[0] >= 0) {
    q = qq;
  } else {
    q = -qq;
  }
  auto const qs = q[0];
  auto const qv = vector_3(q[1], q[2], q[3]);
  auto const qvnorm = norm(qv);
  auto const aanorm = 2.0 * (qvnorm < std::sqrt(0.5) ? std::asin(qvnorm) : std::acos(qs));
  auto const coef = qvnorm < std::sqrt(DBL_EPSILON) ? 2.0 : aanorm / qvnorm;
  auto const aa = coef * qv;
  return aa;
}

// logarithm of a rotation tensor in Special Orthogonal Group(3), as the
// the axis of rotation times the angle of rotation.
OMEGA_H_INLINE Vector<3> axis_angle_from_tensor(Matrix<3, 3> const R) OMEGA_H_NOEXCEPT {
  return axis_angle_from_quaternion(quaternion_from_tensor(R));
}

OMEGA_H_INLINE Vector<1> axis_angle_from_tensor(Matrix<2, 2> const R) OMEGA_H_NOEXCEPT {
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
OMEGA_H_INLINE Vector<4> quaternion_from_axis_angle(Vector<3> const aa) OMEGA_H_NOEXCEPT {
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
  auto const R = 2.0 * outer_product(qv, qv) + 2.0 * qs * cross(qv) + (2.0 * square(qs) - 1.0) * I;
  return R;
}

OMEGA_H_INLINE Matrix<3, 3> tensor_from_axis_angle(Vector<3> const aa) OMEGA_H_NOEXCEPT {
  return tensor_from_quaternion(quaternion_from_axis_angle(aa));
}

OMEGA_H_INLINE Matrix<2, 2> tensor_from_axis_angle(Vector<1> const aa) OMEGA_H_NOEXCEPT {
  auto const theta = aa[0];
  return rotate(theta);
}

OMEGA_H_INLINE Matrix<3, 3> pack_polar(Matrix<3, 3> const spd, Vector<3> const aa) OMEGA_H_NOEXCEPT {
  Matrix<3, 3> packed = spd;
  packed(0,1) = aa[0];
  packed(0,2) = aa[1];
  packed(1,2) = aa[2];
  return packed;
}

OMEGA_H_INLINE Matrix<3, 3> unpack_polar_spd(Matrix<3, 3> const packed) OMEGA_H_NOEXCEPT {
  Matrix<3, 3> spd = packed;
  spd(0,1) = spd(1,0);
  spd(0,2) = spd(2,0);
  spd(1,2) = spd(2,1);
  return spd;
}

OMEGA_H_INLINE Vector<3> unpack_polar_axis_angle(Matrix<3, 3> const packed) OMEGA_H_NOEXCEPT {
  Vector<3> aa;
  aa[0] = packed(0,1);
  aa[1] = packed(0,2);
  aa[2] = packed(1,2);
  return aa;
}

OMEGA_H_INLINE Matrix<2, 2> pack_polar(Matrix<2, 2> const spd, Vector<1> const aa) OMEGA_H_NOEXCEPT {
  Matrix<2, 2> packed = spd;
  packed(0,1) = aa[0];
  return packed;
}

OMEGA_H_INLINE Matrix<2, 2> unpack_polar_spd(Matrix<2, 2> const packed) OMEGA_H_NOEXCEPT {
  Matrix<2, 2> spd = packed;
  spd(0,1) = spd(1,0);
  return spd;
}

OMEGA_H_INLINE Vector<1> unpack_polar_axis_angle(Matrix<2, 2> const packed) OMEGA_H_NOEXCEPT {
  Vector<1> aa;
  aa[0] = packed(0,1);
  return aa;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> log_polar(Matrix<dim, dim> const F) OMEGA_H_NOEXCEPT {
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
OMEGA_H_INLINE Matrix<dim, dim> exp_polar(Matrix<dim, dim> const packed_polar) OMEGA_H_NOEXCEPT {
  auto const r = unpack_polar_axis_angle(packed_polar);
  auto const R = tensor_from_axis_angle(r);
  auto const u = unpack_polar_spd(packed_polar);
  auto const U = exp_spd(u);
  auto const F = R * U;
  return F;
}

OMEGA_H_INLINE Matrix<1, 1> exp_polar(Matrix<1, 1> const log_F) OMEGA_H_NOEXCEPT {
  return matrix_1x1(std::exp(log_F(0, 0)));
}

// this function modifies the axis-angle rotations of a given array of logarithms of polar decompositions
// such that there are no axis-angle vectors which are nearly opposite with angles close to pi.
// this later allows weighted sums of those vectors to give meaningful results and avoids the catastropic
// cancellation case of having (epsilon - pi) and (-epsilon + pi) average to zero.
// tolerance: between 0 and 1, tolerance for opposite pseudo-vectors mapping to rotations close to pi
template <Int dim>
OMEGA_H_INLINE void align_packed_axis_angles(Matrix<dim, dim>* const a, Int const n, Real const tolerance) OMEGA_H_NOEXCEPT {
  auto const alpha = 1.0 - tolerance;
  auto const s_1 = unpack_polar_axis_angle(a[0]);
  auto const s = norm(s_1);
  auto const pi_sq = square(PI);
  auto const pi_2 = 2.0 * PI;
  for (Int i = 1; i < n; ++i) {
    auto s_i = unpack_polar_axis_angle(a[i]);
    if ((s_i * s_1) < (-alpha * pi_sq)) { // pseudo-vectors are nearly opposite with angle close to pi
      s_i = s_i - pi_2 * normalize(s_i);
      a[i] = pack_polar(unpack_polar_spd(a[i]), s_i);
    }
    s = s + norm(s_i);
  }
  if (s > (n * PI)) { // renormalize so that ||s_i|| <= pi
    for (Int i = 0; i < n; ++i) {
      auto s_i = unpack_polar_axis_angle(a[i]);
      s_i = s_i - pi_2 * normalize(s_i);
      a[i] = pack_polar(unpack_polar_spd(a[i]), s_i);
    }
  }
}

OMEGA_H_INLINE void align_packed_axis_angles(Matrix<1, 1>* const, Int const, Real const) OMEGA_H_NOEXCEPT {
}

}  // namespace Omega_h

#endif
