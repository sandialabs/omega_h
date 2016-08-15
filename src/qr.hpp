#ifndef QR_HPP
#define QR_HPP

namespace osh {

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 10.1. Householder QR Factorization

   for k=1 to n
     x = A_{k:m,k}
     v_k = sign(x_1)\|x\|_2 e_1 + x
     v_k = v_k / \|v_k\|_2          <- note this can divide by zero if x={0}
     A_{k:m,k:n} = A_{k:m,k:n} - 2 v_k (v_k^* A_{k:m,k:n}) */

/* the "o" (offset) parameters to householder_vector and reflect_columns
   are there to support hessenberg reduction / tri-diagonalization */

template <Int max_m, Int max_n>
INLINE Vector<max_m> householder_vector(
    Int m, Matrix<max_m, max_n> a, Real anorm, Int k, Int o) {
  Real norm_x = 0;
  for (Int i = k + o; i < m; ++i) norm_x += square(a[k][i]);
  norm_x = sqrt(norm_x);
  /* technically, every matrix has a QR decomposition.
   * if norm_x is close to zero here, the matrix is rank-deficient
   * and we could just skip this reflection and carry forward
   * the rank information.
   * however, all current uses of this code require the matrix
   * to be full-rank, so we can save a bunch of bookkeeping up
   * the stack if we simply assert this here.
   */
  CHECK(norm_x > EPSILON * anorm);
  Vector<max_m> v_k;
  for (Int i = k + o; i < m; ++i) v_k[i] = a[k][i];
  v_k[k + o] += sign(a[k][k + o]) * norm_x;
  Real norm_v_k = 0;
  for (Int i = k + o; i < m; ++i) norm_v_k += square(v_k[i]);
  norm_v_k = sqrt(norm_v_k);
  for (Int i = k + o; i < m; ++i) v_k[i] /= norm_v_k;
  return v_k;
}

template <Int max_m, Int max_n>
INLINE void reflect_columns(
    Int m, Int n, Matrix<max_m, max_n>& a, Vector<max_m> v_k, Int k, Int o) {
  for (Int j = k; j < n; ++j) {
    Real dot = 0;
    for (Int i = k + o; i < m; ++i) dot += a[j][i] * v_k[i];
    for (Int i = k + o; i < m; ++i) a[j][i] -= 2 * dot * v_k[i];
  }
}

template <Int max_m, Int max_n>
struct QRFactorization {
  Few<Vector<max_m>, max_n> v;  // the householder vectors
  Matrix<max_n, max_n> r;
};

template <Int max_m, Int max_n>
INLINE QRFactorization<max_m, max_n> factorize_qr_householder(
    Int m, Int n, Matrix<max_m, max_n> a) {
  Few<Vector<max_m>, max_n> v;
  Real anorm = frobenius_norm(m, n, a);
  for (Int k = 0; k < n; ++k) {
    v[k] = householder_vector(m, a, anorm, k, 0);
    reflect_columns(m, n, a, v[k], k, 0);
  }
  auto r = reduced_r_from_full(n, a);
  return {v, r};
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 10.2. Implicit Calculation of a Product $Q^*b$

   for k=1 to n
     b_{k:m} = b_{k:m} - 2 v_k (v_k^* b_{k:m}) */
template <Int max_m, Int max_n>
INLINE Vector<max_n> implicit_q_trans_b(
    Int m, Int n, Few<Vector<max_m>, max_n> v, Vector<max_m> b) {
  for (Int k = 0; k < n; ++k) {
    Real dot = 0;
    for (Int i = k; i < m; ++i) dot += v[k][i] * b[i];
    for (Int i = k; i < m; ++i) b[i] -= 2 * dot * v[k][i];
  }
  Vector<max_n> qtb;
  for (Int i = 0; i < n; ++i) qtb[i] = b[i];
  return qtb;
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 10.2. Implicit Calculation of a Product $Qx$

   for k=n downto 1
     x_{k:m} = x_{k:m} - 2 v_k (v_k^* b_{k:m}) */
template <Int max_m, Int max_n>
INLINE void implicit_q_x(
    Int m, Int n, Vector<max_m>& x, Few<Vector<max_m>, max_n> v) {
  for (Int k2 = 0; k2 < n; ++k2) {
    Int k = n - k2 - 1;
    Real dot = 0;
    for (Int i = k; i < m; ++i) dot += v[k][i] * x[i];
    for (Int i = k; i < m; ++i) x[i] -= 2 * dot * v[k][i];
  }
}

template <Int max_m, Int max_n>
INLINE Matrix<max_n, max_n> reduced_r_from_full(
    Int n, Matrix<max_m, max_n> fr) {
  Matrix<max_n, max_n> rr;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < n; ++i) rr[j][i] = fr[j][i];
  return rr;
}

template <Int max_m>
INLINE Vector<max_m> solve_upper_triangular(
    Int m, Matrix<max_m, max_m> a, Vector<max_m> b) {
  Vector<max_m> x;
  for (Int ii = 0; ii < m; ++ii) {
    Int i = m - ii - 1;
    x[i] = b[i];
    for (Int j = i + 1; j < m; ++j) x[i] -= a[j][i] * x[j];
    x[i] /= a[i][i];
  }
  return x;
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 11.2 Least Squares via QR factorization

   1. Compute the reduced QR factorization A = \hat{Q}\hat{R}
   2. Compute the vector \hat{Q}^* b
   3. Solve the upper-triangular system \hat{R} x = \hat{Q}^* b for x  */
template <Int max_m, Int max_n>
INLINE Vector<max_n> solve_using_qr(
    Int m, Int n, Matrix<max_m, max_n> a, Vector<max_m> b) {
  auto qr = factorize_qr_householder(m, n, a);
  auto qtb = implicit_q_trans_b(m, n, qr.v, b);
  auto x = solve_upper_triangular(n, qr.r, qtb);
  return x;
}
template <Int max_m, Int max_n>
INLINE Vector<max_n> solve_using_qr(Matrix<max_m, max_n> a, Vector<max_m> b) {
  return solve_using_qr(max_m, max_n, a, b);
}

}  // end namespace osh

#endif
