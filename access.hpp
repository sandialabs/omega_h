template <typename Arr>
INLINE void set_symm(Arr a, Int i, Matrix<3,3> symm) {
  a[i * 6 + 0] = symm[0][0];
  a[i * 6 + 1] = symm[1][1];
  a[i * 6 + 2] = symm[2][2];
  a[i * 6 + 3] = symm[1][0];
  a[i * 6 + 4] = symm[2][1];
  a[i * 6 + 5] = symm[2][0];
}

template <typename Arr>
INLINE Matrix<3,3> get_symm_3(Arr a, Int i) {
  Matrix<3,3> symm;
  symm[0][0] = a[i * 6 + 0];
  symm[1][1] = a[i * 6 + 1];
  symm[2][2] = a[i * 6 + 2];
  symm[1][0] = a[i * 6 + 3];
  symm[2][1] = a[i * 6 + 4];
  symm[2][0] = a[i * 6 + 5];
  symm[0][1] = symm[1][0];
  symm[1][2] = symm[2][1];
  symm[0][2] = symm[2][0];
  return symm;
}

template <Int n, class Arr>
INLINE void set_vec(Arr a, Int i, Vector<n> v) {
  for (Int j = 0; j < n; ++j)
    a[i * n + j] = v[j];
}

template <Int n, class Arr>
INLINE Vector<n> get_vec(Arr a, Int i) {
  Vector<n> v;
  for (Int j = 0; j < n; ++j)
    v[j] = a[i * n + j];
  return v;
}
