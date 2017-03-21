#ifndef OMEGA_H_SCALAR_HPP
#define OMEGA_H_SCALAR_HPP

#include <Omega_h_kokkos.hpp>
#include <cfloat>
#include <climits>
#include <cmath>

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

constexpr OMEGA_H_INLINE Int factorial(Int x) {
  return (x > 1) ? (x * factorial(x - 1)) : 1;
}

OMEGA_H_INLINE Real average(Real a, Real b) { return (a + b) / 2.; }

template <typename T>
constexpr OMEGA_H_INLINE T raise(T x, Int p) {
  return (p > 1) ? (x * raise(x, p - 1)) : x;
}

template <typename T>
constexpr OMEGA_H_INLINE T square(T x) {
  return raise(x, 2);
}

template <typename T>
OMEGA_H_INLINE Real cube(T x) {
  return raise(x, 3);
}

OMEGA_H_INLINE Real sign(Real x) { return (x < 0.0) ? -1.0 : 1.0; }

OMEGA_H_INLINE Real clamp(Real x, Real low, Real high) {
  return min2(max2(x, low), high);
}

template <Int dp>
struct Root;

template <>
struct Root<0> {
  static OMEGA_H_INLINE Real eval(Real) { return 1.0; }
};

template <>
struct Root<1> {
  static OMEGA_H_INLINE Real eval(Real x) { return x; }
};

template <>
struct Root<2> {
  static OMEGA_H_INLINE Real eval(Real x) { return sqrt(x); }
};

template <>
struct Root<3> {
  static OMEGA_H_INLINE Real eval(Real x) { return cbrt(x); }
};

template <Int p>
OMEGA_H_INLINE Real root(Real x) {
  return Root<p>::eval(x);
}

/* compile-time-executable Greatest Common Denominator code */
constexpr OMEGA_H_INLINE Int gcd(Int a, Int b) {
  return (b == 0) ? (a) : (gcd(b, a % b));
}

/* specialization system for raising a Real to a power
 * which is an integer fraction.
 * these templated classes are responsible for reducing
 * the fraction via gdc(), and returning the original value
 * if the fraction is unity.
 */
template <Int np, Int dp, Int cd = gcd(np, dp)>
struct Power : public Power<np / cd, dp / cd> {
  using Power<np / cd, dp / cd>::eval;
  static_assert(cd != 1, "reduced case should be specialized");
};

template <Int np, Int dp>
struct Power<np, dp, 1> {
  static OMEGA_H_INLINE Real eval(Real x) { return root<dp>(raise(x, np)); }
  static_assert(np != dp, "equal case should be specialized");
};

template <Int p>
struct Power<p, p, 1> {
  static OMEGA_H_INLINE Real eval(Real x) { return x; }
};

template <Int np, Int dp>
OMEGA_H_INLINE Real power(Real x) {
  return Power<np, dp>::eval(x);
}

OMEGA_H_INLINE Real power(Real x, Int np, Int dp) {
  switch (np) {
    case 1:
      switch (dp) {
        case 1:
          return power<1, 1>(x);
        case 2:
          return power<1, 2>(x);
        case 3:
          return power<1, 3>(x);
      }
    case 2:
      switch (dp) {
        case 1:
          return power<2, 1>(x);
        case 2:
          return power<2, 2>(x);
        case 3:
          return power<2, 3>(x);
      }
    case 3:
      switch (dp) {
        case 1:
          return power<3, 1>(x);
        case 2:
          return power<3, 2>(x);
        case 3:
          return power<3, 3>(x);
      }
  }
  return -42.0;
}

OMEGA_H_INLINE Real rel_diff_with_floor(Real a, Real b, Real floor = EPSILON) {
  Real am = fabs(a);
  Real bm = fabs(b);
  if (am <= floor && bm <= floor) return 0.0;
  return fabs(b - a) / max2(am, bm);
}

OMEGA_H_INLINE bool are_close(
    Real a, Real b, Real tol = EPSILON, Real floor = EPSILON) {
  return rel_diff_with_floor(a, b, floor) <= tol;
}

}  // namespace Omega_h

#endif
