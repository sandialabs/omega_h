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

template <UInt m, UInt n>
INLINE Vector<m> householder_vector(Matrix<m,n> a, UInt k, UInt o) {
  Real norm_x = 0;
  for (UInt i = k + o; i < m; ++i)
    norm_x += square(a[k][i]);
  norm_x = sqrt(norm_x);
  Vector<m> v_k;
  //if the x vector is nearly zero, use the exact zero vector as the
  //householder vector and avoid extra work and divide by zero below
  if (norm_x < EPSILON) {
    for (UInt i = k + o; i < m; ++i)
      v_k[i] = 0.0;
    return v_k;
  }
  for (UInt i = k + o; i < m; ++i)
    v_k[i] = a[k][i];
  v_k[k + o] += sign(a[k][k + o]) * norm_x;
  Real norm_v_k = 0;
  for (UInt i = k + o; i < m; ++i)
    norm_v_k += square(v_k[i]);
  norm_v_k = sqrt(norm_v_k);
  for (UInt i = k + o; i < m; ++i)
    v_k[i] /= norm_v_k;
  return v_k;
}

template <UInt m, UInt n>
INLINE void reflect_columns(Matrix<m,n>& a, Vector<m> v_k, UInt k, UInt o) {
  for (UInt j = k; j < n; ++j) {
    Real dot = 0;
    for (UInt i = k + o; i < m; ++i)
      dot += a[j][i] * v_k[i];
    for (UInt i = k + o; i < m; ++i)
      a[j][i] -= 2 * dot * v_k[i];
  }
}

template <UInt m, UInt n>
INLINE void factorize_qr_householder(Matrix<m,n>& a,
    Few<Vector<m>, n>& v) {
  for (UInt k = 0; k < n; ++k) {
    v[k] = householder_vector(a, k, 0);
    reflect_columns(a, v[k], k, 0);
  }
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 10.2. Implicit Calculation of a Product $Q^*b$

   for k=1 to n
     b_{k:m} = b_{k:m} - 2 v_k (v_k^* b_{k:m}) */
template <UInt m, UInt n>
INLINE void implicit_q_trans_b(Vector<m>& b,
    Few<Vector<m>, n> v) {
  for (UInt k = 0; k < n; ++k) {
    Real dot = 0;
    for (UInt i = k; i < m; ++i)
      dot += v[k][i] * b[i];
    for (UInt i = k; i < m; ++i)
      b[i] -= 2 * dot * v[k][i];
  }
}

/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 10.2. Implicit Calculation of a Product $Qx$

   for k=n downto 1
     x_{k:m} = x_{k:m} - 2 v_k (v_k^* b_{k:m}) */
template <UInt m, UInt n>
INLINE void implicit_q_x(Vector<m>& x,
    Few<Vector<m>, n> v) {
  for (UInt k2 = 0; k2 < n; ++k2) {
    UInt k = n - k2 - 1;
    Real dot = 0;
    for (UInt i = k; i < m; ++i)
      dot += v[k][i] * x[i];
    for (UInt i = k; i < m; ++i)
      x[i] -= 2 * dot * v[k][i];
  }
}

template <UInt m, UInt n>
INLINE void decompose_qr_reduced(Matrix<m,n> a, Matrix<m,n>& q, Matrix<n,n>& r) {
  Few<Vector<m>, n> v;
  factorize_qr_householder(a, v);
  for (UInt j = 0; j < n; ++j)
    for (UInt i = 0; i < n; ++i)
      r[j][i] = a[j][i];
  q = identity_matrix<m,n>();
  for (UInt j = 0; j < n; ++j)
    implicit_q_x(q[j], v);
}

/* A_{1:m,k+1:m} = A_{1:m,k+1:m} - 2(A_{1:m,k+1:m}v_k)v_k^* */
template <UInt m>
INLINE void reflect_rows(Matrix<m,m>& a, Vector<m> v_k, UInt k) {
  for (UInt i = 0; i < m; ++i) {
    Real dot = 0;
    for (UInt j = k + 1; j < m; ++j)
      dot += a[j][i] * v_k[j];
    for (UInt j = k + 1; j < m; ++j)
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
template <UInt m>
INLINE void householder_hessenberg(Matrix<m,m>& a,
    Few<Vector<m>, m - 2>& v) {
  for (UInt k = 0; k < m - 2; ++k) {
    v[k] = householder_vector(a, k, 1);
    reflect_columns(a, v[k], k, 1);
    reflect_rows(a, v[k], k);
  }
}

template <UInt m>
INLINE
typename std::enable_if<(m > 2)>::type
householder_hessenberg2(Matrix<m,m>& a,
    Matrix<m,m>& q) {
  Few<Vector<m>, m - 2> v;
  householder_hessenberg(a, v);
  q = identity_matrix<m,m>();
  for (UInt j = 0; j < m; ++j)
    implicit_q_x(q[j], v);
}

template <UInt m>
INLINE
typename std::enable_if<!(m > 2)>::type
householder_hessenberg2(Matrix<m,m>&,
    Matrix<m,m>& q) {
  q = identity_matrix<m,m>();
}

template <UInt m>
INLINE bool reduce(Matrix<m,m> a, UInt& n) {
  for (; n >= 2; --n)
    if (fabs(a[n - 2][n - 1]) > EPSILON)
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
template <UInt m>
INLINE Real wilkinson_shift(Matrix<m,m> a, UInt n) {
  auto anm1 = a[n - 2][n - 2];
  auto an   = a[n - 1][n - 1];
  auto bnm1 = a[n - 2][n - 1];
  auto delta  = (anm1 - an) / 2;
  auto denom = fabs(delta) + sqrt(square(delta) + square(bnm1));
  return an - ((sign(delta) * square(bnm1)) / denom);
}

template <UInt m>
INLINE void apply_shift(Matrix<m,m>& a, Real mu) {
  for (UInt i = 0; i < m; ++i)
    a[i][i] -= mu;
}


/* Trefethen, Lloyd N., and David Bau III.
   Numerical linear algebra. Vol. 50. Siam, 1997.
   Algorithm 28.2. "Practical" QR Algorithm

   (Q^{(0)})^T A^{(0)} Q^{(0)}                A^{(0)} is a tridiagonalization of A
   for k = 1,2,...
     Pick a shift \mu^{(k)}                   e.g., Wilkinson shift
     Q^{(k)}R^{(k)} = A^{(k-1)} - \mu^{(k)}I  QR factorization of A^{(k-1)} - \mu^{(k)}I
     A^{(k)} = R^{(k)}Q^{(k)} + \mu^{(k)}I    Recombine factors in reverse order
     If any off-diagonal element A^{(k)}_{j,j+1} is sufficiently close to zero,
       set A_{j,j+1} = A_{j+1,j} = 0 to obtain
       [ A_1   0 ]
       [   0 A_2 ] = A^{(k)}
       and now apply the QR algorithm to A_1 and A_2

  we don't quite implement reduction fully.
  what we do is to check the entry A_{m,m-1} for closeness
  to zero, and in that case try to reduce (A_2 is a scalar).
  even then, we continue to do QR decompositions
  on the full matrix, only the shift computation is affected
  by how far we've reduced */
template <UInt m>
INLINE void qr_eigen(Matrix<m,m>& a, Matrix<m,m>& q,
    UInt max_iters = 100) {
  householder_hessenberg2(a, q);
  UInt n = m;
  for (UInt i = 0; i < max_iters; ++i) {
    if (!reduce(a, n))
      return;
    auto mu = wilkinson_shift(a, n);
    apply_shift(a, mu);
    Matrix<m,m> q_k;
    Matrix<m,m> r_k;
    decompose_qr_reduced(a, q_k, r_k);
    a = r_k * q_k;
    apply_shift(a, -mu);
    q = q * q_k;
  }
  NORETURN();
}

template <UInt m>
INLINE void decompose_eigen_qr(Matrix<m,m> a, Matrix<m,m>& q,
    Matrix<m,m>& l, UInt max_iters = 100) {
  l = a;
  qr_eigen(l, q, max_iters);
}
