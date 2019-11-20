#ifndef OMEGA_H_VECTOR_HPP
#define OMEGA_H_VECTOR_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_few.hpp>

namespace Omega_h {

#ifdef OMEGA_H_USE_KOKKOS

template <Int n>
class Vector : public Few<Real, n> {
 public:
  OMEGA_H_INLINE Vector() {}
  inline Vector(std::initializer_list<Real> l) : Few<Real, n>(l) {}
  OMEGA_H_INLINE void operator=(Vector<n> const& rhs) volatile {
    Few<Real, n>::operator=(rhs);
  }
  OMEGA_H_INLINE void operator=(Vector<n> const& rhs) {
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

#else

template <Int n>
class Vector : public Few<Real, n> {
 public:
  inline Vector() = default;
  inline Vector(std::initializer_list<Real> l) : Few<Real, n>(l) {}
  inline Vector& operator=(Vector const&) = default;
  inline Vector& operator=(Vector&&) = default;
  inline Vector(Vector const&) = default;
  inline Vector(Vector&&) = default;
#define OMEGA_H_VECTOR_AT return Few<Real, n>::operator[](i)
  OMEGA_H_INLINE Real& operator()(Int i) OMEGA_H_NOEXCEPT { OMEGA_H_VECTOR_AT; }
  OMEGA_H_INLINE Real const& operator()(Int i) const OMEGA_H_NOEXCEPT {
    OMEGA_H_VECTOR_AT;
  }
#undef OMEGA_H_VECTOR_AT
};

#endif

template <Int n>
OMEGA_H_INLINE Real* scalar_ptr(Vector<n>& v) {
  return &v[0];
}
template <Int n>
OMEGA_H_INLINE Real const* scalar_ptr(Vector<n> const& v) {
  return &v[0];
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator+(Vector<n> a, Vector<n> b) OMEGA_H_NOEXCEPT {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] + b[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n>& operator+=(
    Vector<n>& a, Vector<n> b) OMEGA_H_NOEXCEPT {
  a = a + b;
  return a;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator-(Vector<n> a, Vector<n> b) OMEGA_H_NOEXCEPT {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] - b[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n>& operator-=(
    Vector<n>& a, Vector<n> b) OMEGA_H_NOEXCEPT {
  a = a - b;
  return a;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator-(Vector<n> a) OMEGA_H_NOEXCEPT {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = -a[i];
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator*(Vector<n> a, Real b)OMEGA_H_NOEXCEPT {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] * b;
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n>& operator*=(Vector<n>& a, Real b) OMEGA_H_NOEXCEPT {
  a = a * b;
  return a;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator*(Real a, Vector<n> b)OMEGA_H_NOEXCEPT {
  return b * a;
}

template <Int n>
OMEGA_H_INLINE Vector<n> operator/(Vector<n> a, Real b) OMEGA_H_NOEXCEPT {
  Vector<n> c;
  for (Int i = 0; i < n; ++i) c[i] = a[i] / b;
  return c;
}

template <Int n>
OMEGA_H_INLINE Vector<n>& operator/=(Vector<n>& a, Real b) OMEGA_H_NOEXCEPT {
  a = a / b;
  return a;
}

template <Int n>
OMEGA_H_INLINE Real operator*(Vector<n> a, Vector<n> b)OMEGA_H_NOEXCEPT {
  return inner_product(a, b);
}

template <Int n>
OMEGA_H_INLINE Real norm_squared(Vector<n> v) OMEGA_H_NOEXCEPT {
  return v * v;
}

template <Int n>
OMEGA_H_INLINE Real norm(Vector<n> v) OMEGA_H_NOEXCEPT {
  return std::sqrt(norm_squared(v));
}

OMEGA_H_INLINE Real norm(Vector<1> v) OMEGA_H_NOEXCEPT {
  return std::abs(v[0]);
}

template <Int n>
OMEGA_H_INLINE Vector<n> normalize(Vector<n> v) OMEGA_H_NOEXCEPT {
  return v / norm(v);
}

OMEGA_H_INLINE Vector<1> vector_1(Real x) OMEGA_H_NOEXCEPT {
  Vector<1> v;
  v[0] = x;
  return v;
}

OMEGA_H_INLINE Vector<2> vector_2(Real x, Real y) OMEGA_H_NOEXCEPT {
  Vector<2> v;
  v[0] = x;
  v[1] = y;
  return v;
}

OMEGA_H_INLINE Vector<3> vector_3(Real x, Real y, Real z) OMEGA_H_NOEXCEPT {
  Vector<3> v;
  v[0] = x;
  v[1] = y;
  v[2] = z;
  return v;
}

template <Int n>
OMEGA_H_INLINE bool are_close(Vector<n> a, Vector<n> b, Real tol = EPSILON,
    Real floor = EPSILON) OMEGA_H_NOEXCEPT {
  for (Int i = 0; i < n; ++i)
    if (!are_close(a[i], b[i], tol, floor)) return false;
  return true;
}

template <Int n>
OMEGA_H_INLINE Vector<n> fill_vector(Real value) OMEGA_H_NOEXCEPT {
  Vector<n> v;
  for (Int i = 0; i < n; ++i) v[i] = value;
  return v;
}

template <Int n>
OMEGA_H_INLINE Vector<n> zero_vector() OMEGA_H_NOEXCEPT {
  return fill_vector<n>(0.0);
}

/* Moore-Penrose pseudo-inverse of a vector */
template <Int n>
OMEGA_H_INLINE Vector<n> pseudo_invert(Vector<n> a) OMEGA_H_NOEXCEPT {
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
OMEGA_H_INLINE Vector<n> positivize(Vector<n> v) OMEGA_H_NOEXCEPT {
  std::uint32_t bits = 0;
  for (Int i = 0; i < n; ++i) bits |= (std::uint32_t(v[i] >= 0.0) << i);
  std::uint32_t neg_bits = (~bits) & ((std::uint32_t(1) << n) - 1);
  if (neg_bits > bits) return v * -1.0;
  return v;
}

OMEGA_H_INLINE Real cross(Vector<2> a, Vector<2> b) OMEGA_H_NOEXCEPT {
  return (a[0] * b[1] - a[1] * b[0]);
}

OMEGA_H_INLINE Vector<3> cross(
    Omega_h::Vector<3> a, Omega_h::Vector<3> b) OMEGA_H_NOEXCEPT {
  return vector_3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]);
}

OMEGA_H_INLINE Vector<2> perp(Vector<2> v) OMEGA_H_NOEXCEPT {
  return vector_2(-v[1], v[0]);
}

template <Int n>
OMEGA_H_INLINE Vector<n> project(Vector<n> a, Vector<n> b) OMEGA_H_NOEXCEPT {
  auto dir = normalize(b);
  auto len = a * dir;
  return len * dir;
}

template <Int n>
OMEGA_H_INLINE Vector<n> reject(Vector<n> a, Vector<n> b) OMEGA_H_NOEXCEPT {
  return a - project(a, b);
}

template <Int n>
OMEGA_H_DEVICE void set_vector(
    Write<Real> const& a, Int i, Vector<n> v) OMEGA_H_NOEXCEPT {
  for (Int j = 0; j < n; ++j) a[i * n + j] = v[j];
}

template <Int n, class Arr>
OMEGA_H_DEVICE Vector<n> get_vector(Arr const& a, Int i) OMEGA_H_NOEXCEPT {
  Vector<n> v;
  for (Int j = 0; j < n; ++j) v[j] = a[i * n + j];
  return v;
}

template <int new_dim, int old_dim>
OMEGA_H_INLINE Vector<new_dim> resize(Vector<old_dim> v) OMEGA_H_NOEXCEPT {
  constexpr int min_dim = Omega_h::min2(new_dim, old_dim);
  Vector<new_dim> v2;
  for (int i = 0; i < min_dim; ++i) v2(i) = v(i);
  for (int i = min_dim; i < new_dim; ++i) v2(i) = 0.0;
  return v2;
}

Reals get_vector_norms(Reals vs, Int dim);
Reals normalize_vectors(Reals vs, Int dim);
Reals dot_vectors(Reals a, Reals b, Int dim);

Reals resize_vectors(Reals vectors, Int old_dim, Int new_dim);

template <Int dim>
Reals repeat_vector(LO n, Vector<dim> v);

extern template Reals repeat_vector(LO n, Vector<1> v);
extern template Reals repeat_vector(LO n, Vector<2> v);
extern template Reals repeat_vector(LO n, Vector<3> v);
extern template Reals repeat_vector(LO n, Vector<4> v);
extern template Reals repeat_vector(LO n, Vector<6> v);
extern template Reals repeat_vector(LO n, Vector<9> v);

}  // namespace Omega_h

#endif
