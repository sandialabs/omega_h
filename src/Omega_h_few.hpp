#ifndef OMEGA_H_FEW_HPP
#define OMEGA_H_FEW_HPP

#include <type_traits>
#include <initializer_list>

#include <Omega_h_kokkos.hpp>

namespace Omega_h {

template <typename T, Int n>
class Few {
  using UninitT = typename std::aligned_storage<sizeof(T), alignof(T)>::type;
  UninitT array_[n];

 public:
  enum { size = n };
  OMEGA_H_INLINE T* data() { return reinterpret_cast<T*>(array_); }
  OMEGA_H_INLINE T const* data() const {
    return reinterpret_cast<T const*>(array_);
  }
  OMEGA_H_INLINE T volatile* data() volatile {
    return reinterpret_cast<T volatile*>(array_);
  }
  OMEGA_H_INLINE T const volatile* data() const volatile {
    return reinterpret_cast<T const volatile*>(array_);
  }
  OMEGA_H_INLINE T& operator[](Int i) { return data()[i]; }
  OMEGA_H_INLINE T const& operator[](Int i) const { return data()[i]; }
  OMEGA_H_INLINE T volatile& operator[](Int i) volatile { return data()[i]; }
  OMEGA_H_INLINE T const volatile& operator[](Int i) const volatile {
    return data()[i];
  }
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (data() + (i++)) T(*it);
    }
  }
  OMEGA_H_INLINE Few() {
    for (Int i = 0; i < n; ++i) new (data() + i) T();
  }
  OMEGA_H_INLINE ~Few() {
    for (Int i = 0; i < n; ++i) (data()[i]).~T();
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) volatile {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE void operator=(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  OMEGA_H_INLINE Few(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
  OMEGA_H_INLINE Few(Few<T, n> const volatile& rhs) {
    for (Int i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
};

}

#endif
