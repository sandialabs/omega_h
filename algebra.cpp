template <UInt m, UInt n>
Matrix<m,n>::Matrix(std::initializer_list<Real> l) {
  UInt k = 0;
  for (Real v : l) {
    (*this)[k % n][k / n] = v;
    ++k;
  }
}
template class Matrix<3,3>;

template <UInt m, UInt n>
std::ostream& operator<<(std::ostream& o, Matrix<m,n> a)
{
  for (UInt i = 0; i < m; ++i) {
    for (UInt j = 0; j < n; ++j)
      o << ' ' << a[j][i];
    o << '\n';
  }
  return o;
}
template std::ostream& operator<<(std::ostream& o, Matrix<3,3> a);
