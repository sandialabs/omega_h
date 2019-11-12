#ifndef OMEGA_H_FEW_HPP
#define OMEGA_H_FEW_HPP

#include <initializer_list>
#include <new>
#include <type_traits>

#include <Omega_h_defines.hpp>
#include <Omega_h_scalar.hpp>

namespace Omega_h {

#ifdef OMEGA_H_USE_KOKKOS

template <typename T, Int n>
class Few {
  T array_[n];

 public:
  using value_type = T;
  OMEGA_H_INLINE T* data() { return array_; }
  OMEGA_H_INLINE T const* data() const { return array_; }
  OMEGA_H_INLINE T volatile* data() volatile { return array_; }
  OMEGA_H_INLINE T const volatile* data() const volatile { return array_; }
  OMEGA_H_INLINE constexpr Int size() const { return n; }
#ifdef OMEGA_H_CHECK_BOUNDS
#define OMEGA_H_FEW_AT                                                         \
  OMEGA_H_CHECK(0 <= i);                                                       \
  OMEGA_H_CHECK(i < size());                                                   \
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
  OMEGA_H_INLINE const T* begin() const OMEGA_H_NOEXCEPT { return array_; }
  OMEGA_H_INLINE const T* end() const OMEGA_H_NOEXCEPT { return array_ + size(); }
  OMEGA_H_INLINE T* begin() OMEGA_H_NOEXCEPT { return array_; }
  OMEGA_H_INLINE T* end() OMEGA_H_NOEXCEPT { return array_ + size(); }
};

#else

template <typename T, Int n>
class Few {
  T array_[n];

 public:
  using value_type = T;
  OMEGA_H_INLINE T* data() OMEGA_H_NOEXCEPT { return array_; }
  OMEGA_H_INLINE T const* data() const OMEGA_H_NOEXCEPT { return array_; }
  OMEGA_H_INLINE constexpr Int size() const { return n; }
#ifdef OMEGA_H_CHECK_BOUNDS
#define OMEGA_H_FEW_AT                                                         \
  OMEGA_H_CHECK(0 <= i);                                                       \
  OMEGA_H_CHECK(i < size());                                                   \
  return array_[i]
#else
#define OMEGA_H_FEW_AT return array_[i]
#endif
  OMEGA_H_INLINE T& operator[](Int i) OMEGA_H_NOEXCEPT { OMEGA_H_FEW_AT; }
  OMEGA_H_INLINE T const& operator[](Int i) const OMEGA_H_NOEXCEPT {
    OMEGA_H_FEW_AT;
  }
#undef OMEGA_H_FEW_AT
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (array_ + (i++)) T(*it);
    }
  }
  inline Few() = default;
  inline ~Few() = default;
  inline Few(Few<T, n> const& rhs) = default;
  inline Few(Few<T, n>&& rhs) = default;
  inline Few& operator=(Few const& rhs) = default;
  inline Few& operator=(Few&& rhs) = default;
  OMEGA_H_INLINE const T* begin() const OMEGA_H_NOEXCEPT { return array_; }
  OMEGA_H_INLINE const T* end() const OMEGA_H_NOEXCEPT { return array_ + size(); }
  OMEGA_H_INLINE T* begin() OMEGA_H_NOEXCEPT { return array_; }
  OMEGA_H_INLINE T* end() OMEGA_H_NOEXCEPT { return array_ + size(); }
};

#endif

template <Int capacity, typename T>
OMEGA_H_INLINE void add_unique(
    Few<T, capacity>& stack, Int& n, T e) OMEGA_H_NOEXCEPT {
  for (Int i = 0; i < n; ++i)
    if (stack[i] == e) return;
  stack[n++] = e;
}

template <Int n, typename T>
OMEGA_H_INLINE T average(Few<T, n> x) OMEGA_H_NOEXCEPT {
  auto avg = x[0];
  for (Int i = 1; i < n; ++i) avg = avg + x[i];
  return avg / n;
}

template <Int n, typename T, typename Op>
OMEGA_H_INLINE T reduce(Few<T, n> x, Op op) OMEGA_H_NOEXCEPT {
  auto out = x[0];
  for (Int i = 1; i < n; ++i) out = op(out, x[i]);
  return out;
}

#if !(defined(OMEGA_H_USE_CUDA) && defined(__clang__))
template <Int n, typename T>
OMEGA_H_INLINE decltype(std::declval<T>() * std::declval<T>()) inner_product(
    Few<T, n> a, Few<T, n> b) OMEGA_H_NOEXCEPT {
  auto out = a[0] * b[0];
  for (Int i = 1; i < n; ++i) out = out + (a[i] * b[i]);
  return out;
}
#else
template <Int n, typename T>
OMEGA_H_INLINE decltype(T() * T()) inner_product(
    Few<T, n> a, Few<T, n> b) OMEGA_H_NOEXCEPT {
  auto out = a[0] * b[0];
  for (Int i = 1; i < n; ++i) out = out + (a[i] * b[i]);
  return out;
}
#endif

}  // namespace Omega_h

#endif
