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

template <Int m, Int n>
OSH_INLINE bool householder_vector(Matrix<m,n> a, Real anorm,
    Int k, Int o,
    Vector<m>& v_k) {
  Real norm_x = 0;
  for (Int i = k + o; i < m; ++i)
    norm_x += square(a[k][i]);
  norm_x = sqrt(norm_x);
  //if the x vector is nearly zero, use the exact zero vector as the
  //householder vector and avoid extra work and divide by zero below
  if (norm_x <= EPSILON * anorm) {
    for (Int i = k + o; i < m; ++i)
      v_k[i] = 0.0;
    return false;
  }
  for (Int i = k + o; i < m; ++i)
    v_k[i] = a[k][i];
  v_k[k + o] += sign(a[k][k + o]) * norm_x;
  Real norm_v_k = 0;
  for (Int i = k + o; i < m; ++i)
    norm_v_k += square(v_k[i]);
  norm_v_k = sqrt(norm_v_k);
  for (Int i = k + o; i < m; ++i)
    v_k[i] /= norm_v_k;
  return true;
}

template <Int m, Int n>
OSH_INLINE void reflect_columns(Matrix<m,n>& a, Vector<m> v_k, Int k, Int o) {
  for (Int j = k; j < n; ++j) {
    Real dot = 0;
    for (Int i = k + o; i < m; ++i)
      dot += a[j][i] * v_k[i];
    for (Int i = k + o; i < m; ++i)
      a[j][i] -= 2 * dot * v_k[i];
  }
}

template <Int m, Int n>
OSH_INLINE Int factorize_qr_householder(Matrix<m,n>& a,
    Few<Vector<m>, n>& v) {
  Real anorm = frobenius_norm(a);
  Int rank = 0;
  for (Int k = 0; k < n; ++k) {
    rank += householder_vector(a, anorm, k, 0, v[k]);
    reflect_columns(a, v[k], k, 0);
  }
  return rank;
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 10.2. Implicit Calculation of a Product $Q^*b$

   for k=1 to n
     b_{k:m} = b_{k:m} - 2 v_k (v_k^* b_{k:m}) */
template <Int m, Int n>
OSH_INLINE void implicit_q_trans_b(Vector<m>& b,
    Few<Vector<m>, n> v) {
  for (Int k = 0; k < n; ++k) {
    Real dot = 0;
    for (Int i = k; i < m; ++i)
      dot += v[k][i] * b[i];
    for (Int i = k; i < m; ++i)
      b[i] -= 2 * dot * v[k][i];
  }
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 10.2. Implicit Calculation of a Product $Qx$

   for k=n downto 1
     x_{k:m} = x_{k:m} - 2 v_k (v_k^* b_{k:m}) */
template <Int m, Int n>
OSH_INLINE void implicit_q_x(Vector<m>& x,
    Few<Vector<m>, n> v) {
  for (Int k2 = 0; k2 < n; ++k2) {
    Int k = n - k2 - 1;
    Real dot = 0;
    for (Int i = k; i < m; ++i)
      dot += v[k][i] * x[i];
    for (Int i = k; i < m; ++i)
      x[i] -= 2 * dot * v[k][i];
  }
}

template <Int m, Int n>
OSH_INLINE Matrix<n,n> reduced_r_from_full(Matrix<m,n> fr) {
  Matrix<n,n> rr;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < n; ++i)
      rr[j][i] = fr[j][i];
  return rr;
}

template <Int m, Int n>
OSH_INLINE Int decompose_qr_reduced(Matrix<m,n> a, Matrix<m,n>& q, Matrix<n,n>& r) {
  Few<Vector<m>, n> v;
  Int rank = factorize_qr_householder(a, v);
  r = reduced_r_from_full(a);
  q = identity_matrix<m,n>();
  for (Int j = 0; j < n; ++j)
    implicit_q_x(q[j], v);
  return rank;
}

/* A_{1:m,k+1:m} = A_{1:m,k+1:m} - 2(A_{1:m,k+1:m}v_k)v_k^* */
template <Int m>
OSH_INLINE void reflect_rows(Matrix<m,m>& a, Vector<m> v_k, Int k) {
  for (Int i = 0; i < m; ++i) {
    Real dot = 0;
    for (Int j = k + 1; j < m; ++j)
      dot += a[j][i] * v_k[j];
    for (Int j = k + 1; j < m; ++j)
      a[j][i] -= 2 * dot * v_k[j];
  }
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 26.1. Householder Reduction to Hessenberg Form

   (note the similarity to Algorithm 10.1, hence the code reuse)

   for k=1 to m - 2
     x = A_{k+1:m,k}
     v_k = sign(x_1)\|x\|_2 e_1 + x
     v_k = v_k / \|v_k\|_2
     A_{k+1:m,k:m} = A_{k+1:m,k:m} - 2 v_k (v_k^* A_{k+1:m,k:m})
     A_{1:m,k+1:m} = A_{1:m,k+1:m} - 2(A_{1:m,k+1:m}v_k)v_k^* */
template <Int m>
OSH_INLINE void householder_hessenberg(Matrix<m,m>& a,
    Few<Vector<m>, m - 2>& v) {
  Real anorm = frobenius_norm(a);
  for (Int k = 0; k < m - 2; ++k) {
    householder_vector(a, anorm, k, 1, v[k]);
    reflect_columns(a, v[k], k, 1);
    reflect_rows(a, v[k], k);
  }
}

template <Int m>
OSH_INLINE
typename std::enable_if<(m > 2)>::type
householder_hessenberg2(Matrix<m,m>& a,
    Matrix<m,m>& q) {
  Few<Vector<m>, m - 2> v;
  householder_hessenberg(a, v);
  q = identity_matrix<m,m>();
  for (Int j = 0; j < m; ++j)
    implicit_q_x(q[j], v);
}

template <Int m>
OSH_INLINE
typename std::enable_if<!(m > 2)>::type
householder_hessenberg2(Matrix<m,m>&,
    Matrix<m,m>& q) {
  q = identity_matrix<m,m>();
}

template <Int m>
OSH_INLINE bool reduce(Matrix<m,m> a, Real anorm, Int& n) {
  for (; n >= 2; --n)
    if (fabs(a[n - 2][n - 1]) > EPSILON * anorm)
      return true;
  return false;
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Equation 29.8. Wilkinson Shift

   let B denote the lower-rightmost 2x2 submatrix of A^{(k)}

   B = [ a_{m-1} b_{m-1} ]
       [ b_{m-1} a_{m}   ]

   \mu = (a_m - sign(\delta)b_{m-1}^2) /
         (|\delta| + \sqrt{\delta^2 + b_{m-1}^2})

   where \delta = (a_{m-1} - a_m) / 2    */
template <Int m>
OSH_INLINE Real wilkinson_shift(Matrix<m,m> a, Int n) {
  auto anm1 = a[n - 2][n - 2];
  auto an   = a[n - 1][n - 1];
  auto bnm1 = a[n - 2][n - 1];
  auto delta  = (anm1 - an) / 2;
  auto denom = fabs(delta) + sqrt(square(delta) + square(bnm1));
  return an - ((sign(delta) * square(bnm1)) / denom);
}

template <Int m>
OSH_INLINE void apply_shift(Matrix<m,m>& a, Real mu) {
  subtract_from_diag(a, mu);
}

template <Int m>
OSH_INLINE Vector<m> solve_upper_triangular(Matrix<m,m> a, Vector<m> b) {
  Vector<m> x;
  for (Int ii = 0; ii < m; ++ii) {
    Int i = m - ii - 1;
    x[i] = b[i];
    for (Int j = i + 1; j < m; ++j)
      x[i] -= a[j][i] * x[j];
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
template <Int m, Int n>
OSH_INLINE bool solve_least_squares_qr(Matrix<m,n> a, Vector<m> b,
    Vector<n>& x) {
  Few<Vector<m>, n> v;
  Int rank = factorize_qr_householder(a, v);
  if (rank != n)
    return false;
  Matrix<n,n> r = reduced_r_from_full(a);
  Vector<m> qtb_full = b;
  implicit_q_trans_b(qtb_full, v);
  Vector<n> qtb;
  for (Int i = 0; i < n; ++i)
    qtb[i] = qtb_full[i];
  x = solve_upper_triangular(r, qtb);
  return true;
}
