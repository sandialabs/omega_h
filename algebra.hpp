INLINE Real average(Real a, Real b) { return (a + b) / 2.; }

INLINE Real square(Real x) { return x * x; }

INLINE Real cube(Real x) { return x * x * x; }

INLINE Real sign(Real x) { return (x < 0.0) ? -1.0 : 1.0; }

INLINE bool are_close(Real a, Real b,
    Real tol = EPSILON, Real floor = EPSILON) {
  Real am = fabs(a);
  Real bm = fabs(b);
  if (am < floor && bm < floor)
    return true;
  if ((fabs(b - a) / max2(am, bm)) <= tol)
    return true;
  return false;
}

template <Int n>
class Vector : public Few<Real, n> {
  public:
    INLINE Vector() {}
    Vector(std::initializer_list<Real> l);
    INLINE void operator=(Vector<n> const& rhs) volatile {
      Few<Real, n>::operator=(rhs);
    }
    INLINE Vector(Vector<n> const& rhs):
      Few<Real, n>(rhs) {
    }
    INLINE Vector(const volatile Vector<n>& rhs):
      Few<Real, n>(rhs) {
    }
};

template <Int n>
Vector<n>::Vector(std::initializer_list<Real> l) {
  Int i = 0;
  for (Real v : l)
    (*this)[i++] = v;
}

template <Int n>
std::ostream& operator<<(std::ostream& o, Vector<n> v) {
  for (Int i = 0; i < n; ++i)
    o << ' ' << v[i];
  o << '\n';
  return o;
}

template <Int n>
INLINE Vector<n> operator*(Vector<n> a, Real b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i)
    c[i] = a[i] * b;
  return c;
}

template <Int n>
INLINE Vector<n> operator/(Vector<n> a, Real b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i)
    c[i] = a[i] / b;
  return c;
}

template <Int n>
INLINE Real operator*(Vector<n> a, Vector<n> b) {
  Real c = a[0] * b[0];
  for (Int i = 1; i < n; ++i)
    c += a[i] * b[i];
  return c;
}

template <Int n>
INLINE Real norm_squared(Vector<n> v) {
  return v * v;
}

template <Int n>
INLINE Real norm(Vector<n> v) {
  return sqrt(norm_squared(v));
}

template <Int n>
INLINE Vector<n> normalize(Vector<n> v) {
  return v / norm(v);
}

template <Int n>
INLINE Vector<n> operator+(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i)
    c[i] = a[i] + b[i];
  return c;
}

template <Int n>
INLINE Vector<n> operator-(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i)
    c[i] = a[i] - b[i];
  return c;
}

template <Int n>
INLINE bool are_close(Vector<n> a, Vector<n> b,
    Real tol = EPSILON, Real floor = EPSILON) {
  for (Int i = 0; i < n; ++i)
    if (!are_close(a[i], b[i], tol, floor))
      return false;
  return true;
}

/* column-first storage and indexing !!! */
template <Int m, Int n>
class Matrix : public Few<Vector<m>, n> {
  public:
    INLINE Matrix() {}
    Matrix(std::initializer_list<Real> l);
    INLINE void operator=(Matrix<m,n> const& rhs) volatile {
      Few<Vector<m>, n>::operator=(rhs);
    }
    INLINE Matrix(Matrix<m,n> const& rhs):
      Few<Vector<m>, n>(rhs) {
    }
    INLINE Matrix(const volatile Matrix<m,n>& rhs):
      Few<Vector<m>, n>(rhs) {
    }
};

template <Int m, Int n>
Matrix<m,n>::Matrix(std::initializer_list<Real> l) {
  Int k = 0;
  for (Real v : l) {
    (*this)[k % n][k / n] = v;
    ++k;
  }
}

template <Int m, Int n>
std::ostream& operator<<(std::ostream& o, Matrix<m,n> a)
{
  for (Int i = 0; i < m; ++i) {
    for (Int j = 0; j < n; ++j)
      o << ' ' << a[j][i];
    o << '\n';
  }
  return o;
}

template <Int m, Int n>
INLINE Matrix<m,n> identity_matrix() {
  Matrix<m,n> a;
  for (Int j = 0; j < n; ++j)
  for (Int i = 0; i < m; ++i)
    a[j][i] = (i == j);
  return a;
}

template <Int m, Int n>
INLINE Vector<m> operator*(Matrix<m,n> a, Vector<n> b) {
  Vector<m> c = a[0] * b[0];
  for (Int j = 1; j < n; ++j)
    c = c + a[j] * b[j];
  return c;
}

template <Int m, Int n, Int p>
INLINE Matrix<m,n> operator*(Matrix<m,p> a, Matrix<p,n> b) {
  Matrix<m,n> c;
  for (Int j = 0; j < n; ++j)
    c[j] = a * b[j];
  return c;
}

template <Int m, Int n>
INLINE Matrix<m,n> operator*(Matrix<m,n> a, Real b) {
  Matrix<m,n> c;
  for (Int j = 0; j < n; ++j)
    c[j] = a[j] * b;
  return c;
}

template <Int m, Int n>
INLINE Matrix<m,n> operator*(Real a, Matrix<m,n> b) {
  return b * a;
}

template <Int m, Int n>
INLINE Matrix<m,n> operator/(Matrix<m,n> a, Real b) {
  Matrix<m,n> c;
  for (Int j = 0; j < n; ++j)
    c[j] = a[j] / b;
  return c;
}

template <Int m, Int n>
INLINE Matrix<m,n> operator+(Matrix<m,n> a, Matrix<m,n> b) {
  Matrix<m,n> c;
  for (Int j = 0; j < n; ++j)
    c[j] = a[j] + b[j];
  return c;
}

template <Int m, Int n>
INLINE Matrix<m,n> operator-(Matrix<m,n> a, Matrix<m,n> b) {
  Matrix<m,n> c;
  for (Int j = 0; j < n; ++j)
    c[j] = a[j] - b[j];
  return c;
}

template <Int m, Int n>
INLINE Matrix<n,m> transpose(Matrix<m,n> a) {
  Matrix<n,m> b;
  for (Int i = 0; i < m; ++i)
  for (Int j = 0; j < n; ++j)
    b[i][j] = a[j][i];
  return b;
}

template <Int m, Int n>
INLINE Real max_norm(Matrix<m,n> a) {
  Real x = 0.0;
  for (Int j = 0; j < n; ++j)
  for (Int i = 0; i < m; ++i)
    x = max2(x, fabs(a[j][i]));
  return x;
}

template <Int m, Int n>
INLINE Real frobenius_norm(Matrix<m,n> a) {
  Real x = 0.0;
  for (Int j = 0; j < n; ++j)
  for (Int i = 0; i < m; ++i)
    x += square(a[j][i]);
  return sqrt(x);
}

template <Int m, Int n>
INLINE bool are_close(Matrix<m,n> a, Matrix<m,n> b,
    Real tol = EPSILON, Real floor = EPSILON) {
  for (Int j = 0; j < n; ++j)
    if (!are_close(a[j], b[j], tol, floor))
      return false;
  return true;
}

template <Int m, Int n>
INLINE Matrix<m,n> tensor_product(Vector<m> a, Vector<n> b) {
  Matrix<m,n> c;
  for (Int j = 0; j < n; ++j)
  for (Int i = 0; i < m; ++i)
    c[j][i] = a[i] * b[j];
  return c;
}

template <Int m>
INLINE Real trace(Matrix<m,m> a) {
  Real t = a[0][0];
  for (Int i = 1; i < m; ++i)
    t += a[i][i];
  return t;
}

template <Int m>
INLINE Vector<m> diagonal(Matrix<m,m> a) {
  Vector<m> v;
  for (Int i = 0; i < m; ++i)
    v[i] = a[i][i];
  return v;
}

template <Int m>
INLINE Matrix<m,m> diagonal(Vector<m> v) {
  Matrix<m,m> a;
  for (Int i = 0; i < m; ++i)
  for (Int j = i + 1; j < m; ++j)
    a[i][j] = a[j][i] = 0.0;
  for (Int i = 0; i < m; ++i)
    a[i][i] = v[i];
  return a;
}

template <Int n>
INLINE Vector<n> zero_vector() {
  Vector<n> v;
  for (Int i = 0; i < n; ++i)
    v[i] = 0.0;
  return v;
}

template <Int m, Int n>
INLINE Matrix<m,n> zero_matrix() {
  Matrix<m,n> a;
  for (Int j = 0; j < n; ++j)
    a[j] = zero_vector<m>();
  return a;
}

template <Int m>
INLINE void subtract_from_diag(Matrix<m,m>& a, Real mu) {
  for (Int i = 0; i < m; ++i)
    a[i][i] -= mu;
}
