#ifndef OMEGA_H_MATH_HPP
#define OMEGA_H_MATH_HPP

#include <cfloat>
#include <climits>
#include <cmath>

#include "Omega_h.hpp"

namespace Omega_h {

template <typename T>
struct ArithTraits;

template <>
struct ArithTraits<unsigned char> {
  static OMEGA_H_INLINE unsigned char max() { return UCHAR_MAX; }
  static OMEGA_H_INLINE unsigned char min() { return 0; }
};

template <>
struct ArithTraits<signed char> {
  static OMEGA_H_INLINE signed char max() { return SCHAR_MAX; }
  static OMEGA_H_INLINE signed char min() { return SCHAR_MIN; }
};

template <>
struct ArithTraits<unsigned int> {
  static OMEGA_H_INLINE unsigned int max() { return UINT_MAX; }
  static OMEGA_H_INLINE unsigned int min() { return 0; }
};

template <>
struct ArithTraits<int> {
  static OMEGA_H_INLINE int max() { return INT_MAX; }
  static OMEGA_H_INLINE int min() { return INT_MIN; }
};

template <>
struct ArithTraits<unsigned long> {
  static OMEGA_H_INLINE unsigned long max() { return ULONG_MAX; }
  static OMEGA_H_INLINE unsigned long min() { return 0; }
};

template <>
struct ArithTraits<signed long> {
  static OMEGA_H_INLINE signed long max() { return LONG_MAX; }
  static OMEGA_H_INLINE signed long min() { return LONG_MIN; }
};

template <>
struct ArithTraits<unsigned long long> {
  static OMEGA_H_INLINE unsigned long long max() { return ULLONG_MAX; }
  static OMEGA_H_INLINE unsigned long long min() { return 0; }
};

template <>
struct ArithTraits<signed long long> {
  static OMEGA_H_INLINE signed long long max() { return LLONG_MAX; }
  static OMEGA_H_INLINE signed long long min() { return LLONG_MIN; }
};

template <>
struct ArithTraits<double> {
  static OMEGA_H_INLINE double max() { return DBL_MAX; }
  static OMEGA_H_INLINE double min() { return -DBL_MAX; }
};

template <typename T>
OMEGA_H_INLINE T square(T x) {
  return x * x;
}

template <Int n>
class Vector : public Few<Real, n> {
 public:
  OMEGA_H_INLINE Vector() {}
  inline Vector(std::initializer_list<Real> l) : Few<Real, n>(l) {}
  OMEGA_H_INLINE void operator=(Vector<n> const& rhs) volatile {
    Few<Real, n>::operator=(rhs);
  }
  OMEGA_H_INLINE Vector(Vector<n> const& rhs) : Few<Real, n>(rhs) {}
  OMEGA_H_INLINE Vector(const volatile Vector<n>& rhs) : Few<Real, n>(rhs) {}
};

template <Int n>
OMEGA_H_INLINE Vector<n> operator+(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] + b[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator-(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] - b[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator*(Vector<n> a, Real b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] * b;
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator*(Real a, Vector<n> b) {
  return b * a;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator/(Vector<n> a, Real b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] / b;
  return c;
}

template <Int n>
OMEGA_H_INLINE Real operator*(Vector<n> a, Vector<n> b) {
  Real c = a[0] * b[0];
  for (Int i = 1; i < n; ++i) c += a[i] * b[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Real norm_squared(Vector<n> v) {
  return v * v;
}

template <Int n>
OMEGA_H_INLINE Real norm(Vector<n> v) {
  return sqrt(norm_squared(v));
}

template <Int n>
OMEGA_H_INLINE Vector<n> normalize(Vector<n> v) {
  return v / norm(v);
}

OMEGA_H_INLINE Vector<2> vector_2(Real x, Real y) {
  Vector<2> v;
  v[0] = x;
  v[1] = y;
  return v;
}

OMEGA_H_INLINE Vector<3> vector_3(Real x, Real y, Real z) {
  Vector<3> v;
  v[0] = x;
  v[1] = y;
  v[2] = z;
  return v;
}

/* column-first storage and indexing !!! */
template <Int m, Int n>
class Matrix : public Few<Vector<m>, n> {
 public:
  OMEGA_H_INLINE Matrix() {}
  /* these constructors accept the matrix in
   * row-first order for convenience.
   */
  inline Matrix(std::initializer_list<Vector<m>> l) : Few<Vector<m>, n>(l) {}
  inline Matrix(std::initializer_list<Real> l);
  OMEGA_H_INLINE void operator=(Matrix<m, n> const& rhs) volatile {
    Few<Vector<m>, n>::operator=(rhs);
  }
  OMEGA_H_INLINE Matrix(Matrix<m, n> const& rhs) : Few<Vector<m>, n>(rhs) {}
  OMEGA_H_INLINE Matrix(const volatile Matrix<m, n>& rhs)
      : Few<Vector<m>, n>(rhs) {}
};

template <Int m, Int n>
OMEGA_H_INLINE Vector<m> operator*(Matrix<m, n> a, Vector<n> b) {
  Vector<m> c = a[0] * b[0];
  for (Int j = 1; j < n; ++j) c = c + a[j] * b[j];
  return c;
}

template <Int m, Int n, Int p>
OMEGA_H_INLINE Matrix<m, n> operator*(Matrix<m, p> a, Matrix<p, n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a * b[j];
  return c;
}

template <Int m, Int n>
OMEGA_H_INLINE Matrix<n, m> transpose(Matrix<m, n> a) {
  Matrix<n, m> b;
  for (Int i = 0; i < m; ++i)
    for (Int j = 0; j < n; ++j) b[i][j] = a[j][i];
  return b;
}

template <Int max_m, Int max_n>
OMEGA_H_INLINE Matrix<max_m, max_n> identity_matrix(Int m, Int n) {
  Matrix<max_m, max_n> a;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) a[j][i] = (i == j);
  return a;
}

template <Int max_m, Int max_n>
OMEGA_H_INLINE Matrix<max_m, max_n> identity_matrix() {
  return identity_matrix<max_m, max_n>(max_m, max_n);
}

OMEGA_H_INLINE Matrix<3, 3> matrix_3x3(
    Real a, Real b, Real c, Real d, Real e, Real f, Real g, Real h, Real i) {
  Matrix<3, 3> o;
  o[0] = vector_3(a, d, g);
  o[1] = vector_3(b, e, h);
  o[2] = vector_3(c, f, i);
  return o;
}

OMEGA_H_INLINE Matrix<3, 3> cross(Vector<3> a) {
  return matrix_3x3(0, -a[2], a[1], a[2], 0, -a[0], -a[1], a[0], 0);
}

OMEGA_H_INLINE Vector<3> cross(Vector<3> a, Vector<3> b) {
  return vector_3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]);
}

template <Int m>
OMEGA_H_INLINE Matrix<m, m> diagonal(Vector<m> v) {
  Matrix<m, m> a;
  for (Int i = 0; i < m; ++i)
    for (Int j = i + 1; j < m; ++j) a[i][j] = a[j][i] = 0.0;
  for (Int i = 0; i < m; ++i) a[i][i] = v[i];
  return a;
}

template <Int n>
OMEGA_H_DEVICE void set_vector(Write<Real> const& a, Int i, Vector<n> v) {
  for (Int j = 0; j < n; ++j) a[i * n + j] = v[j];
}

OMEGA_H_INLINE constexpr Int symm_dofs(Int dim) {
  return ((dim + 1) * dim) / 2;
}

OMEGA_H_INLINE Vector<3> symm2vector(Matrix<2, 2> symm) {
  Vector<3> v;
  v[0] = symm[0][0];
  v[1] = symm[1][1];
  v[2] = symm[1][0];
  return v;
}

OMEGA_H_INLINE Vector<6> symm2vector(Matrix<3, 3> symm) {
  Vector<6> v;
  v[0] = symm[0][0];
  v[1] = symm[1][1];
  v[2] = symm[2][2];
  v[3] = symm[1][0];
  v[4] = symm[2][1];
  v[5] = symm[2][0];
  return v;
}

template <Int n>
OMEGA_H_DEVICE void set_symm(Write<Real> const& a, Int i, Matrix<n, n> symm) {
  set_vector(a, i, symm2vector(symm));
}

template <Int dim>
OMEGA_H_INLINE Vector<dim> metric_eigenvalues(Vector<dim> h) {
  Vector<dim> l;
  for (Int i = 0; i < dim; ++i) l[i] = 1.0 / square(h[i]);
  return l;
}

template <Int dim>
OMEGA_H_INLINE Matrix<dim, dim> compose_metric(
    Matrix<dim, dim> r, Vector<dim> h) {
  auto l = metric_eigenvalues(h);
  return r * diagonal(l) * transpose(r);
}

template <Int dim>
Reals repeat_symm(LO n, Matrix<dim, dim> symm);
extern template Reals repeat_symm(LO n, Matrix<3, 3> symm);
extern template Reals repeat_symm(LO n, Matrix<2, 2> symm);

}  // end namespace Omega_h

#endif
