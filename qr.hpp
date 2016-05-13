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
