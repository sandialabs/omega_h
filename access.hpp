template <Int m, class Arr>
INLINE void set_symm(Arr a, Int i, Matrix<m,m> symm) {
  a[i * 6 + 0] = symm[0][0];
  a[i * 6 + 1] = symm[1][1];
  a[i * 6 + 2] = symm[2][2];
  a[i * 6 + 3] = symm[1][0];
  a[i * 6 + 4] = symm[2][1];
  a[i * 6 + 5] = symm[2][0];
}

template <Int m, class Arr>
INLINE Matrix<m,m> get_symm(Arr a, Int i) {
  Matrix<m,m> symm;
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
