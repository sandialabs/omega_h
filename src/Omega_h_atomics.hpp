#ifndef OMEGA_H_ATOMICS_HPP
#define OMEGA_H_ATOMICS_HPP

namespace Omega_h {

// TODO: remove the whole disable system once the Kokkos February 2018 version
// is released

template <bool Enable>
struct Atomics;

#ifdef OMEGA_H_USE_KOKKOSCORE
template <>
struct Atomics<true> {
  template <typename T>
  static OMEGA_H_INLINE void increment(volatile T* const dest) {
    Kokkos::atomic_increment(dest);
  }
  template <typename T>
  static OMEGA_H_INLINE void add(volatile T* const dest, const T val) {
    Kokkos::atomic_add(dest, val);
  }
  template <typename T>
  static OMEGA_H_INLINE T fetch_add(volatile T* const dest, const T val) {
    return Kokkos::atomic_fetch_add(dest, val);
  }
};
#endif

template <>
struct Atomics<false> {
  template <typename T>
  static OMEGA_H_INLINE void increment(volatile T* const dest) {
    ++(*dest);
  }
  template <typename T>
  static OMEGA_H_INLINE void add(volatile T* const dest, const T val) {
    *dest += val;
  }
  template <typename T>
  static OMEGA_H_INLINE T fetch_add(volatile T* const dest, const T val) {
    T tmp = *dest;
    *dest += val;
    return tmp;
  }
};

#ifdef OMEGA_H_USE_KOKKOSCORE
constexpr bool enable_atomics =
    !std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Serial>::value;
#else
constexpr bool enable_atomics = false;
#endif

template <typename T>
OMEGA_H_INLINE void atomic_increment(volatile T* const dest) {
  Atomics<enable_atomics>::increment<T>(dest);
}

template <typename T>
OMEGA_H_INLINE void atomic_add(volatile T* const dest, const T val) {
  Atomics<enable_atomics>::add<T>(dest, val);
}

template <typename T>
OMEGA_H_INLINE T atomic_fetch_add(volatile T* const dest, const T val) {
  return Atomics<enable_atomics>::fetch_add<T>(dest, val);
}

}  // end namespace Omega_h

#endif
