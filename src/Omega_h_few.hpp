#ifndef OMEGA_H_FEW_HPP
#define OMEGA_H_FEW_HPP

#include <initializer_list>
#include <new>
#include <type_traits>

#include <Omega_h_defines.hpp>
#include <Omega_h_scalar.hpp>

namespace Omega_h {

template <typename T, Int n>
class Few {
  /* Can't use this because volatile forces us to make the class look non-POD
   * static_assert(std::is_pod<T>::value, "T must be POD for Few<T>");
   */
  T array_[n];

 public:
  using value_type = T;
  enum { size = n };
  OMEGA_H_INLINE T* data() { return array_; }
  OMEGA_H_INLINE T const* data() const { return array_; }
  OMEGA_H_INLINE T volatile* data() volatile { return array_; }
  OMEGA_H_INLINE T const volatile* data() const volatile { return array_; }
#ifdef OMEGA_H_CHECK_BOUNDS
#define OMEGA_H_FEW_AT                                                         \
  OMEGA_H_CHECK(0 <= i);                                                       \
  OMEGA_H_CHECK(i < size);                                                     \
  return array_[i]
#else
#define OMEGA_H_FEW_AT return array_[i]
#endif
  OMEGA_H_INLINE T& operator[](Int i) { OMEGA_H_FEW_AT; }
  OMEGA_H_INLINE T const& operator[](Int i) const { OMEGA_H_FEW_AT; }
  OMEGA_H_INLINE T volatile& operator[](Int i) volatile { OMEGA_H_FEW_AT; }
  OMEGA_H_INLINE T const volatile& operator[](Int i) const volatile {
    OMEGA_H_FEW_AT;
  }
#undef OMEGA_H_FEW_AT
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (array_ + (i++)) T(*it);
    }
  }
  OMEGA_H_INLINE Few() {}
  OMEGA_H_INLINE ~Few() {}
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) volatile {
    for (Int i = 0; i < n; ++i) array_[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) array_[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) array_[i] = rhs[i];
  }
  OMEGA_H_INLINE Few(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) new (array_ + i) T(rhs[i]);
  }
  OMEGA_H_INLINE Few(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) new (array_ + i) T(rhs[i]);
  }
};

template <Int capacity, typename T>
OMEGA_H_INLINE void add_unique(Few<T, capacity>& stack, Int& n, T e) {
  for (Int i = 0; i < n; ++i)
    if (stack[i] == e) return;
  stack[n++] = e;
}

template <Int n, typename T>
OMEGA_H_INLINE T average(Few<T, n> x) {
  auto avg = x[0];
  for (Int i = 1; i < n; ++i) avg = avg + x[i];
  return avg / n;
}

template <Int n, typename T>
OMEGA_H_INLINE T minimum(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = min2(out, x[i]);
  return out;
}

template <Int n, typename T>
OMEGA_H_INLINE T maximum(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = max2(out, x[i]);
  return out;
}

template <Int n_max, typename T>
OMEGA_H_INLINE T maximum_magnitude(Few<T, n_max> x, Int n) {
  auto out = x[0];
  auto max_mag = std::abs(x[0]);
  for (Int i = 1; i < n; ++i) {
    auto mag = std::abs(x[i]);
    if (mag > max_mag) {
      max_mag = mag;
      out = x[i];
    }
  }
  return out;
}

template <Int n, typename T>
OMEGA_H_INLINE T sum(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = out + x[i];
  return out;
}

template <Int n, typename T>
OMEGA_H_INLINE T product(Few<T, n> x) {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = out * x[i];
  return out;
}

#if !(defined(OMEGA_H_USE_CUDA) && defined(__clang__))
template <Int n, typename T>
OMEGA_H_INLINE decltype(std::declval<T>() * std::declval<T>()) inner_product(
    Few<T, n> a, Few<T, n> b) {
  auto out = a[0] * b[0];
  for (Int i = 1; i < n; ++i) out = out + (a[i] * b[i]);
  return out;
}
#else
template <Int n, typename T>
OMEGA_H_INLINE decltype(T() * T()) inner_product(Few<T, n> a, Few<T, n> b) {
  auto out = a[0] * b[0];
  for (Int i = 1; i < n; ++i) out = out + (a[i] * b[i]);
  return out;
}
#endif

}  // namespace Omega_h

#endif
