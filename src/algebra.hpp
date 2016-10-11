#ifndef ALGEBRA_HPP
#define ALGEBRA_HPP

#include <cmath>

#include "Omega_h_math.hpp"
#include "internal.hpp"

namespace Omega_h {

INLINE Real average(Real a, Real b) { return (a + b) / 2.; }

INLINE Real cube(Real x) { return x * x * x; }

INLINE Real sign(Real x) { return (x < 0.0) ? -1.0 : 1.0; }

template <Int p>
INLINE Real raise(Real x) {
  Real out = x;
  for (Int i = 1; i < p; ++i) out *= x;
  return out;
}

template <Int dp>
struct Root;

template <>
struct Root<0> {
  static INLINE Real eval(Real) { return 1.0; }
};

template <>
struct Root<1> {
  static INLINE Real eval(Real x) { return x; }
};

template <>
struct Root<2> {
  static INLINE Real eval(Real x) { return sqrt(x); }
};

template <>
struct Root<3> {
  static INLINE Real eval(Real x) { return cbrt(x); }
};

template <Int np, Int dp>
struct Power {
  static INLINE Real eval(Real x) { return Root<dp>::eval(raise<np>(x)); }
  static_assert(np != dp, "equal case should be specialized");
};

template <Int p>
struct Power<p, p> {
  static INLINE Real eval(Real x) { return x; }
};

template <Int np, Int dp>
INLINE Real power(Real x) {
  return Power<np, dp>::eval(x);
}

INLINE Real rel_diff_with_floor(Real a, Real b, Real floor = EPSILON) {
  Real am = fabs(a);
  Real bm = fabs(b);
  if (am <= floor && bm <= floor) return 0.0;
  return fabs(b - a) / max2(am, bm);
}

INLINE bool are_close(
    Real a, Real b, Real tol = EPSILON, Real floor = EPSILON) {
  return rel_diff_with_floor(a, b, floor) <= tol;
}

template <Int n, typename T>
INLINE T average(Few<T, n> x) {
  auto avg = x[0];
  for (Int i = 1; i < n; ++i) avg = avg + x[i];
  return avg / n;
}

template <Int n>
INLINE Vector<n> operator-(Vector<n> a) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = -a[i];
  return c;
}

template <Int n>
INLINE bool are_close(
    Vector<n> a, Vector<n> b, Real tol = EPSILON, Real floor = EPSILON) {
  for (Int i = 0; i < n; ++i)
    if (!are_close(a[i], b[i], tol, floor)) return false;
  return true;
}

template <Int n>
INLINE Vector<n> zero_vector() {
  Vector<n> v;
  for (Int i = 0; i < n; ++i) v[i] = 0.0;
  return v;
}

/* Moore-Penrose pseudo-inverse of a vector */
template <Int n>
INLINE Vector<n> pseudo_invert(Vector<n> a) {
  auto nsq = a * a;
  if (nsq < EPSILON) return zero_vector<n>();
  return a / nsq;
}

template <Int m, Int n>
INLINE Matrix<m, n> operator*(Matrix<m, n> a, Real b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] * b;
  return c;
}

template <Int m, Int n>
INLINE Matrix<m, n> operator*(Real a, Matrix<m, n> b) {
  return b * a;
}

template <Int m, Int n>
INLINE Matrix<m, n> operator/(Matrix<m, n> a, Real b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] / b;
  return c;
}

template <Int m, Int n>
INLINE Matrix<m, n> operator+(Matrix<m, n> a, Matrix<m, n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] + b[j];
  return c;
}

template <Int m, Int n>
INLINE Matrix<m, n> operator-(Matrix<m, n> a, Matrix<m, n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j) c[j] = a[j] - b[j];
  return c;
}

template <Int m, Int n>
INLINE Real max_norm(Matrix<m, n> a) {
  Real x = 0.0;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) x = max2(x, fabs(a[j][i]));
  return x;
}

template <Int max_m, Int max_n>
INLINE Real frobenius_norm(Int m, Int n, Matrix<max_m, max_n> a) {
  Real x = 0.0;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) x += square(a[j][i]);
  return sqrt(x);
}

template <Int m, Int n>
INLINE bool are_close(
    Matrix<m, n> a, Matrix<m, n> b, Real tol = EPSILON, Real floor = EPSILON) {
  for (Int j = 0; j < n; ++j)
    if (!are_close(a[j], b[j], tol, floor)) return false;
  return true;
}

template <Int m, Int n>
INLINE Matrix<m, n> tensor_product(Vector<m> a, Vector<n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) c[j][i] = a[i] * b[j];
  return c;
}

template <Int m>
INLINE Real trace(Matrix<m, m> a) {
  Real t = a[0][0];
  for (Int i = 1; i < m; ++i) t += a[i][i];
  return t;
}

template <Int m>
INLINE Vector<m> diagonal(Matrix<m, m> a) {
  Vector<m> v;
  for (Int i = 0; i < m; ++i) v[i] = a[i][i];
  return v;
}

template <Int m, Int n>
INLINE Matrix<m, n> zero_matrix() {
  Matrix<m, n> a;
  for (Int j = 0; j < n; ++j) a[j] = zero_vector<m>();
  return a;
}

template <Int m>
INLINE void subtract_from_diag(Matrix<m, m>& a, Real mu) {
  for (Int i = 0; i < m; ++i) a[i][i] -= mu;
}

/* a function to disambiguate a unit vector
   from its negative. we treat the signs of
   the components as bits of an integer,
   and negate the components if the resulting
   bit pattern makes a larger integer */
template <Int n>
INLINE Vector<n> positivize(Vector<n> v) {
  std::uint32_t bits = 0;
  for (Int i = 0; i < n; ++i) bits |= (std::uint32_t(v[i] >= 0.0) << i);
  std::uint32_t neg_bits = (~bits) & ((std::uint32_t(1) << n) - 1);
  if (neg_bits > bits) return v * -1.0;
  return v;
}

Reals normalize_vectors(Reals vs, Int dim);
Reals interpolate_between(Reals a, Reals b, Real t);

}  // end namespace Omega_h

#endif
