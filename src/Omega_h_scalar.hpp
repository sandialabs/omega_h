#ifndef OMEGA_H_SCALAR_HPP
#define OMEGA_H_SCALAR_HPP

#include <climits>
#include <cfloat>
#include <Omega_h_kokkos.hpp>

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
constexpr OMEGA_H_INLINE T square(T x) {
  return x * x;
}

template <typename T>
OMEGA_H_INLINE T max2(T a, T b) {
  return (b > a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE T min2(T a, T b) {
  return (b < a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE void swap2(T& a, T& b) {
  T c = a;
  a = b;
  b = c;
}

}

#endif
