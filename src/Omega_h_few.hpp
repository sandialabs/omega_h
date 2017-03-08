#ifndef OMEGA_H_FEW_HPP
#define OMEGA_H_FEW_HPP

#include <type_traits>
#include <initializer_list>

#include <Omega_h_kokkos.hpp>
#include <Omega_h_defines.hpp>

namespace Omega_h {

template <typename T, Int n>
class Few {
/* Can't use this because volatile forces us to make the class look non-POD
 * static_assert(std::is_pod<T>::value, "T must be POD for Few<T>");
 */
  T array_[n];

 public:
  enum { size = n };
  OMEGA_H_INLINE T const* data() const {
    return array_;
  }
  OMEGA_H_INLINE T volatile* data() volatile {
    return array_;
  }
  OMEGA_H_INLINE T const volatile* data() const volatile {
    return array_;
  }
  OMEGA_H_INLINE T& operator[](Int i) { return array_[i]; }
  OMEGA_H_INLINE T const& operator[](Int i) const { return array_[i]; }
  OMEGA_H_INLINE T volatile& operator[](Int i) volatile { return array_[i]; }
  OMEGA_H_INLINE T const volatile& operator[](Int i) const volatile {
    return array_[i];
  }
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (array_ + (i++)) T(*it);
    }
  }
  OMEGA_H_INLINE Few() {
    for (Int i = 0; i < n; ++i) array_[i] = T();
  }
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
    for (Int i = 0; i < n; ++i) array_[i] = rhs[i];
  }
  OMEGA_H_INLINE Few(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) array_[i] = rhs[i];
  }
};

}

#endif
