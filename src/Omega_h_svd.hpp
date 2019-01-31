#ifndef OMEGA_H_SVD_HPP
#define OMEGA_H_SVD_HPP

#include <Omega_h_eigen.hpp>

namespace Omega_h {

//
// Givens rotation. [c, -s; s, c] [a; b] = [r; 0]
// \param a, b
// \return c, s
//

struct Givens {
  Real c;
  Real s;
};

OMEGA_H_INLINE
Givens givens(Real const a, Real const b) noexcept {
  auto c = 1.0;
  auto s = 0.0;
  if (b != 0.0) {
    if (std::abs(b) > std::abs(a)) {
      auto const t = -a / b;
      s = 1.0 / std::sqrt(1.0 + t * t);
      c = t * s;
    } else {
      auto const t = -b / a;
      c = 1.0 / std::sqrt(1.0 + t * t);
      s = t * c;
    }
  }
  return {c, s};
}

template <Int dim>
struct SVD {
  Matrix<dim, dim> U;
  Matrix<dim, dim> S;
  Matrix<dim, dim> V;
};

//
// Singular value decomposition (SVD) for 2x2
// bidiagonal matrix. Used for general 2x2 SVD.
// Adapted from LAPAPCK's DLASV2, Netlib's dlasv2.c
// and LBNL computational crystallography toolbox
// \param f, g, h where A = [f, g; 0, h]
// \return \f$ A = USV^T\f$
//
OMEGA_H_INLINE
SVD<2> svd_bidiagonal(Real f, Real const g, Real h) noexcept {
  auto fa = std::abs(f);
  auto ga = std::abs(g);
  auto ha = std::abs(h);
  auto s0 = 0.0;
  auto s1 = 0.0;
  auto cu = 1.0;
  auto su = 0.0;
  auto cv = 1.0;
  auto sv = 0.0;
  auto const swap_diag = (ha > fa);
  if (swap_diag == true) {
    swap2(fa, ha);
    swap2(f, h);
  }
  // diagonal matrix
  if (ga == 0.0) {
    s1 = ha;
    s0 = fa;
  } else if (ga > fa && fa / ga < DBL_EPSILON) {
    // case of very large ga
    s0 = ga;
    s1 = ha > 1.0 ? Real(fa / (ga / ha)) : Real((fa / ga) * ha);
    cu = 1.0;
    su = h / g;
    cv = f / g;
    sv = 1.0;
  } else {
    // normal case
    auto const d = fa - ha;
    auto const l = d / fa;   // l \in [0,1]
    auto const m = g / f;    // m \in (-1/macheps, 1/macheps)
    auto const t = 2.0 - l;  // t \in [1,2]
    auto const mm = m * m;
    auto const tt = t * t;
    auto const s = std::sqrt(tt + mm);  // s \in [1,1 + 1/macheps]
    auto const r = ((l != 0.0) ? (std::sqrt(l * l + mm))
                               : (std::abs(m)));  // r \in [0,1 + 1/macheps]
    auto const a = 0.5 * (s + r);                 // a \in [1,1 + |m|]
    s1 = ha / a;
    s0 = fa * a;
    // Compute singular vectors
    Real tau;  // second assignment to T in DLASV2
    if (mm != 0.0) {
      tau = (m / (s + t) + m / (r + l)) * (1.0 + a);
    } else {
      // note that m is very tiny
      tau = (l == 0.0) ? (std::copysign(2.0, f) * std::copysign(1.0, g))
                       : (g / std::copysign(d, f) + m / t);
    }
    auto const lv =
        std::sqrt(tau * tau + 4.0);  // second assignment to L in DLASV2
    cv = 2.0 / lv;
    sv = tau / lv;
    cu = (cv + sv * m) / a;
    su = (h / f) * sv / a;
  }
  // Fix signs of singular values in accordance to sign of singular vectors
  s0 = std::copysign(s0, f);
  s1 = std::copysign(s1, h);
  if (swap_diag == true) {
    swap2(cu, sv);
    swap2(su, cv);
  }
  auto const U = matrix_2x2(cu, -su, su, cu);
  auto const S = matrix_2x2(s0, 0.0, 0.0, s1);
  auto const V = matrix_2x2(cv, -sv, sv, cv);
  return {U, S, V};
}

OMEGA_H_INLINE
SVD<2> svd_2x2(Matrix<2, 2> const A) OMEGA_H_NOEXCEPT {
  // First compute a givens rotation to eliminate 1,0 entry in tensor
  auto const cs = givens(A(0, 0), A(1, 0));
  auto const c = cs.c;
  auto const s = cs.s;
  auto const R = matrix_2x2(c, -s, s, c);
  auto const B = R * A;
  // B is bidiagonal. Use specialized algorithm to compute its SVD
  auto const XSV = svd_bidiagonal(B(0, 0), B(0, 1), B(1, 1));
  auto const X = XSV.U;
  auto const S = XSV.S;
  auto const V = XSV.V;
  // Complete general 2x2 SVD with givens rotation calculated above
  auto const U = transpose(R) * X;
  return {U, S, V};
}

//
// R^N singular value decomposition (SVD)
// \param A tensor
// \return \f$ A = USV^T\f$
//
template <Int dim>
OMEGA_H_INLINE SVD<dim> decompose_svd(
    Matrix<dim, dim> const A) OMEGA_H_NOEXCEPT {
  // Scale first
  auto const norm_a = norm(A);
  auto const scale = norm_a > 0.0 ? norm_a : Real(1.0);
  auto S = A / scale;
  auto U = identity_matrix<dim, dim>();
  auto V = identity_matrix<dim, dim>();
  auto off = norm_off_diag(S);
  auto const tol = DBL_EPSILON;
  Int const max_iter = 2048;
  Int num_iter = 0;
  while (off > tol && num_iter < max_iter) {
    // Find largest off-diagonal entry
    auto const pq = arg_max_off_diag(S);
    auto p = pq[0];
    auto q = pq[1];
    if (p > q) swap2(p, q);
    // Obtain left and right Givens rotations by using 2x2 SVD
    auto const Spq = matrix_2x2(S(p, p), S(p, q), S(q, p), S(q, q));
    auto const LDR = svd_2x2(Spq);
    auto const L = LDR.U;
    auto const R = LDR.V;
    auto const cl = L(0, 0);
    auto const sl = L(0, 1);
    auto const cr = R(0, 0);
    auto const sr = (sign(R(0, 1)) == sign(R(1, 0))) ? (-R(0, 1)) : (R(0, 1));
    // Apply both Givens rotations to matrices
    // that are converging to singular values and singular vectors
    S = givens_left(cl, sl, p, q, S);
    S = givens_right(cr, sr, p, q, S);
    U = givens_right(cl, sl, p, q, U);
    V = givens_left(cr, sr, p, q, V);
    off = norm_off_diag(S);
    ++num_iter;
  }
  // Fix signs for entries in the diagonal matrix S
  // that are negative
  for (Int i = 0; i < dim; ++i) {
    if (S(i, i) < 0.0) {
      S(i, i) = -S(i, i);
      for (Int j = 0; j < dim; ++j) {
        U(j, i) = -U(j, i);
      }
    }
  }
  S *= scale;
  return {U, S, V};
}

}  // namespace Omega_h

#endif
