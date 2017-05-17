#ifndef OMEGA_H_VECTOR_HPP
#define OMEGA_H_VECTOR_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_few.hpp>

namespace Omega_h {

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
#define OMEGA_H_VECTOR_AT return Few<Real, n>::operator[](i)
  OMEGA_H_INLINE Real& operator()(Int i) { OMEGA_H_VECTOR_AT; }
  OMEGA_H_INLINE Real const& operator()(Int i) const { OMEGA_H_VECTOR_AT; }
  OMEGA_H_INLINE Real volatile& operator()(Int i) volatile {
    OMEGA_H_VECTOR_AT;
  }
  OMEGA_H_INLINE Real const volatile& operator()(Int i) const volatile {
    OMEGA_H_VECTOR_AT;
  }
#undef OMEGA_H_VECTOR_AT
};

template <Int n>
OMEGA_H_INLINE Real* scalar_ptr(Vector<n>& v) {
  return &v[0];
}
template <Int n>
OMEGA_H_INLINE Real const* scalar_ptr(Vector<n> const& v) {
  return &v[0];
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator+(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] + b[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n>& operator+=(Vector<n>& a, Vector<n> b) {
  a = a + b;
  return a;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator-(Vector<n> a, Vector<n> b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] - b[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n>& operator-=(Vector<n>& a, Vector<n> b) {
  a = a - b;
  return a;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator-(Vector<n> a) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = -a[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator*(Vector<n> a, Real b) {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] * b;
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n>& operator*=(Vector<n>& a, Real b) {
  a = a * b;
  return a;
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
OMEGA_H_INLINE Vector<n>& operator/=(Vector<n>& a, Real b) {
  a = a / b;
  return a;
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

OMEGA_H_INLINE Vector<1> vector_1(Real x) {
  Vector<1> v;
  v[0] = x;
  return v;
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

template <Int n>
OMEGA_H_INLINE bool are_close(
    Vector<n> a, Vector<n> b, Real tol = EPSILON, Real floor = EPSILON) {
  for (Int i = 0; i < n; ++i)
    if (!are_close(a[i], b[i], tol, floor)) return false;
  return true;
}

template <Int n>
OMEGA_H_INLINE Vector<n> fill_vector(Real value) {
  Vector<n> v;
  for (Int i = 0; i < n; ++i) v[i] = value;
  return v;
}

template <Int n>
OMEGA_H_INLINE Vector<n> zero_vector() {
  return fill_vector<n>(0.0);
}

/* Moore-Penrose pseudo-inverse of a vector */
template <Int n>
OMEGA_H_INLINE Vector<n> pseudo_invert(Vector<n> a) {
  auto nsq = a * a;
  if (nsq < EPSILON) return zero_vector<n>();
  return a / nsq;
}

/* a function to disambiguate a unit vector
   from its negative. we treat the signs of
   the components as bits of an integer,
   and negate the components if the resulting
   bit pattern makes a larger integer */
template <Int n>
OMEGA_H_INLINE Vector<n> positivize(Vector<n> v) {
  std::uint32_t bits = 0;
  for (Int i = 0; i < n; ++i) bits |= (std::uint32_t(v[i] >= 0.0) << i);
  std::uint32_t neg_bits = (~bits) & ((std::uint32_t(1) << n) - 1);
  if (neg_bits > bits) return v * -1.0;
  return v;
}

OMEGA_H_INLINE Real cross(Vector<2> a, Vector<2> b) {
  return (a[0] * b[1] - a[1] * b[0]);
}

OMEGA_H_INLINE Vector<3> cross(Omega_h::Vector<3> a, Omega_h::Vector<3> b) {
  return vector_3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]);
}

OMEGA_H_INLINE Vector<2> perp(Vector<2> v) { return vector_2(-v[1], v[0]); }

template <Int n>
OMEGA_H_DEVICE void set_vector(Write<Real> const& a, Int i, Vector<n> v) {
  for (Int j = 0; j < n; ++j) a[i * n + j] = v[j];
}

template <Int n, class Arr>
OMEGA_H_DEVICE Vector<n> get_vector(Arr const& a, Int i) {
  Vector<n> v;
  for (Int j = 0; j < n; ++j) v[j] = a[i * n + j];
  return v;
}

Reals get_vector_norms(Reals vs, Int dim);
Reals normalize_vectors(Reals vs, Int dim);

Reals resize_vectors(Reals vectors, Int old_dim, Int new_dim);

}  // namespace Omega_h

#endif
