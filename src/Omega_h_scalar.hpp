#ifndef OMEGA_H_SCALAR_HPP
#define OMEGA_H_SCALAR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_fail.hpp>
#include <cfloat>
#include <climits>
#include <cmath>
#include <utility>

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
constexpr OMEGA_H_INLINE T max2(T a, T b) {
  return (a < b) ? (b) : (a);
}

template <typename T>
constexpr OMEGA_H_INLINE T min2(T a, T b) {
  return (b < a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE void swap2(T& a, T& b) {
  T c(std::move(a));
  a = std::move(b);
  b = std::move(c);
}

template <typename T>
constexpr OMEGA_H_INLINE T factorial(T x) {
  return (x > 1) ? (x * factorial(x - 1)) : 1;
}

template <typename T>
constexpr OMEGA_H_INLINE T average(T a, T b) {
  return (a + b) / 2;
}

template <Int p, typename T>
struct Raise {
  static_assert(p >= 0, "negative power not allowed in Raise!");
  static constexpr OMEGA_H_INLINE T eval(T x) {
    return x * Raise<p - 1, T>::eval(x);
  }
};

template <typename T>
struct Raise<0, T> {
  static constexpr OMEGA_H_INLINE T eval(T) { return 1; }
};

template <Int p, typename T>
constexpr OMEGA_H_INLINE T raise(T x) {
  return Raise<p, T>::eval(x);
}

template <typename T>
constexpr OMEGA_H_INLINE T square(T x) {
  return raise<2, T>(x);
}

template <typename T>
OMEGA_H_INLINE T cube(T x) {
  return raise<3, T>(x);
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
  static OMEGA_H_INLINE Real eval(Real x) { return std::sqrt(x); }
};

template <>
struct Root<3> {
  static OMEGA_H_INLINE Real eval(Real x) { return std::cbrt(x); }
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
  static OMEGA_H_INLINE Real eval(Real x) { return root<dp>(raise<np>(x)); }
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
  Real am = std::abs(a);
  Real bm = std::abs(b);
  if (am <= floor && bm <= floor) return 0.0;
  return std::abs(b - a) / max2(am, bm);
}

OMEGA_H_INLINE bool are_close(
    Real a, Real b, Real tol = EPSILON, Real floor = EPSILON) {
  return rel_diff_with_floor(a, b, floor) <= tol;
}

template <typename T>
T divide_no_remainder(T a, T b) {
  OMEGA_H_CHECK(b != 0);
  OMEGA_H_CHECK(a % b == 0);
  return a / b;
}

}  // namespace Omega_h

#endif
