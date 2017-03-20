#ifndef ALGEBRA_HPP
#define ALGEBRA_HPP

#include <cmath>

#include "Omega_h_math.hpp"
#include "internal.hpp"

namespace Omega_h {

template <Int n, typename T>
INLINE T average(Few<T, n> x) {
  auto avg = x[0];
  for (Int i = 1; i < n; ++i) avg = avg + x[i];
  return avg / n;
}

template <Int n, typename T>
INLINE T minimum(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = min2(out, x[i]);
  return out;
}

template <Int n, typename T>
INLINE T maximum(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = max2(out, x[i]);
  return out;
}

template <Int n, typename T>
INLINE T sum(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = out + x[i];
  return out;
}

template <Int n, typename T>
INLINE T product(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = out * x[i];
  return out;
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
INLINE Matrix<m, n> outer_product(Vector<m> a, Vector<n> b) {
  Matrix<m, n> c;
  for (Int j = 0; j < n; ++j)
    for (Int i = 0; i < m; ++i) c[j][i] = a[i] * b[j];
  return c;
}

template <Int m>
INLINE Real trace(Matrix<m, m> a) __attribute__((pure));
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

INLINE Real determinant(Matrix<1, 1> m) {
  return m[0][0];
}

INLINE Real determinant(Matrix<2, 2> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return a * d - b * c;
}

INLINE Real determinant(Matrix<3, 3> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[2][0];
  Real d = m[0][1];
  Real e = m[1][1];
  Real f = m[2][1];
  Real g = m[0][2];
  Real h = m[1][2];
  Real i = m[2][2];
  return (a * e * i) + (b * f * g) + (c * d * h) - (c * e * g) - (b * d * i) -
         (a * f * h);
}

INLINE Matrix<1, 1> invert(Matrix<1, 1> m) { return matrix_1x1(1.0 / m[0][0]); }

INLINE Matrix<2, 2> invert(Matrix<2, 2> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return matrix_2x2(d, -b, -c, a) / determinant(m);
}

INLINE Matrix<3, 3> invert(Matrix<3, 3> a) {
  Matrix<3, 3> b;
  b[0] = cross(a[1], a[2]);
  b[1] = cross(a[2], a[0]);
  b[2] = cross(a[0], a[1]);
  return transpose(b) / determinant(a);
}

/* Moore-Penrose pseudo-inverse of a vector */
template <Int n>
INLINE Vector<n> pseudo_invert(Vector<n> a) {
  auto nsq = a * a;
  if (nsq < EPSILON) return zero_vector<n>();
  return a / nsq;
}

template <Int m, Int n>
INLINE typename std::enable_if<(n < m), Matrix<n, m>>::type pseudo_invert(
    Matrix<m, n> a) {
  auto at = transpose(a);
  return invert(at * a) * at;
}

template <Int m, Int n>
INLINE typename std::enable_if<(n > m), Matrix<n, m>>::type pseudo_invert(
    Matrix<m, n> a) {
  auto at = transpose(a);
  return at * invert(a * at);
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

Reals get_vector_norms(Reals vs, Int dim);
Reals normalize_vectors(Reals vs, Int dim);
Reals interpolate_between(Reals a, Reals b, Real t);

}  // end namespace Omega_h

#endif
