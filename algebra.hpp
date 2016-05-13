INLINE Real square(Real x) { return x * x; }
INLINE Real sign(Real x) { return (x < 0.0) ? -1.0 : 1.0; }
INLINE bool are_close(Real a, Real b, Real tol, Real floor) {
  Real am = fabs(a);
  Real bm = fabs(b);
  if (am < floor && bm < floor)
    return true;
  if (fabs(b - a) / max2(am, bm) <= tol)
    return true;
  return false;
}

template <UInt n>
class Vector : public Few<Real, n> {
};

template <UInt n>
std::ostream& operator<<(std::ostream& o, Vector<n> v)
{
  for (UInt i = 0; i < n; ++i)
    o << ' ' << v[i];
  o << '\n';
  return o;
}

template <UInt n>
Vector<n> operator*(Vector<n> a, Real b)
{
  Vector<n> c;
  for (UInt i = 0; i < n; ++i)
    c[i] = a[i] * b;
  return c;
}

template <UInt n>
Vector<n> operator+(Vector<n> a, Vector<n> b)
{
  Vector<n> c;
  for (UInt i = 0; i < n; ++i)
    c[i] = a[i] + b[i];
  return c;
}

template <UInt n>
INLINE bool are_close(Vector<n> a, Vector<n> b) {
  for (UInt i = 0; i < n; ++i)
    if (!are_close(a[i], b[i]))
      return false;
  return true;
}

/* column-first storage and indexing !!! */
template <UInt m, UInt n>
class Matrix : public Few<Vector<m>, n> {
  public:
    Matrix() {}
    Matrix(std::initializer_list<Real> l) {
      UInt k = 0;
      for (Real v : l) {
        (*this)[k % n][k / n] = v;
        ++k;
      }
    }
};

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

template <UInt m, UInt n>
Vector<m> operator*(Matrix<m,n> a, Vector<n> b)
{
  Vector<m> c = a[0] * b[0];
  for (UInt j = 1; j < n; ++j)
    c = c + a[j] * b[j];
  return c;
}

template <UInt m, UInt n, UInt p>
Matrix<m,n> operator*(Matrix<m,p> a, Matrix<p,n> b)
{
  Matrix<m,n> c;
  for (UInt j = 0; j < n; ++j)
    c[j] = a * b[j];
  return c;
}

template <UInt m, UInt n>
INLINE bool are_close(Matrix<m,n> a, Matrix<m,n> b) {
  for (UInt j = 0; j < n; ++j)
    if (!are_close(a[j], b[j]))
      return false;
  return true;
}
