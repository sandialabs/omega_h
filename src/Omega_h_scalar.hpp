#ifndef OMEGA_H_SCALAR_HPP
#define OMEGA_H_SCALAR_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_fail.hpp>
#include <cfloat>
#include <climits>
#include <cmath>
#include <utility>

namespace Omega_h {

/* certain operations, including Kokkos reductions and writing
   to std::cout, don't behave as desired on std::int8_t.
   This class is just responsible for raising std::int8_t to std::int32_t */
template <typename T>
struct Promoted {
  typedef T type;
};

template <>
struct Promoted<I8> {
  typedef I32 type;
};

template <typename T>
using promoted_t = typename Promoted<T>::type;

template <typename T>
struct ArithTraits;

template <>
struct ArithTraits<unsigned char> {
  static OMEGA_H_INLINE unsigned char max() noexcept { return UCHAR_MAX; }
  static OMEGA_H_INLINE unsigned char min() noexcept { return 0; }
};

template <>
struct ArithTraits<signed char> {
  static constexpr OMEGA_H_INLINE signed char max() noexcept {
    return SCHAR_MAX;
  }
  static constexpr OMEGA_H_INLINE signed char min() noexcept {
    return SCHAR_MIN;
  }
};

template <>
struct ArithTraits<unsigned int> {
  static constexpr OMEGA_H_INLINE unsigned int max() noexcept {
    return UINT_MAX;
  }
  static constexpr OMEGA_H_INLINE unsigned int min() noexcept { return 0; }
};

template <>
struct ArithTraits<int> {
  static constexpr OMEGA_H_INLINE int max() noexcept { return INT_MAX; }
  static constexpr OMEGA_H_INLINE int min() noexcept { return INT_MIN; }
};

template <>
struct ArithTraits<unsigned long> {
  static constexpr OMEGA_H_INLINE unsigned long max() noexcept {
    return ULONG_MAX;
  }
  static constexpr OMEGA_H_INLINE unsigned long min() noexcept { return 0; }
};

template <>
struct ArithTraits<signed long> {
  static constexpr OMEGA_H_INLINE signed long max() noexcept {
    return LONG_MAX;
  }
  static constexpr OMEGA_H_INLINE signed long min() noexcept {
    return LONG_MIN;
  }
};

template <>
struct ArithTraits<unsigned long long> {
  static constexpr OMEGA_H_INLINE unsigned long long max() noexcept {
    return ULLONG_MAX;
  }
  static constexpr OMEGA_H_INLINE unsigned long long min() noexcept {
    return 0;
  }
};

template <>
struct ArithTraits<signed long long> {
  static constexpr OMEGA_H_INLINE signed long long max() noexcept {
    return LLONG_MAX;
  }
  static constexpr OMEGA_H_INLINE signed long long min() noexcept {
    return LLONG_MIN;
  }
};

template <>
struct ArithTraits<double> {
  static constexpr OMEGA_H_INLINE double max() noexcept { return DBL_MAX; }
  static constexpr OMEGA_H_INLINE double min() noexcept { return -DBL_MAX; }
};

template <typename T>
constexpr OMEGA_H_INLINE T max2(T a, T b) noexcept {
  return (a < b) ? (b) : (a);
}

template <typename T>
constexpr OMEGA_H_INLINE T min2(T a, T b) noexcept {
  return (b < a) ? (b) : (a);
}

template <typename T>
OMEGA_H_INLINE void swap2(T& a, T& b) noexcept {
  T c(std::move(a));
  a = std::move(b);
  b = std::move(c);
}

template <typename T>
constexpr OMEGA_H_INLINE_BIG T factorial(T x) noexcept {
  return (x > 1) ? (x * factorial(x - 1)) : 1;
}

template <typename T>
constexpr OMEGA_H_INLINE T average(T a, T b) noexcept {
  return (a + b) / 2;
}

template <Int p, typename T>
struct Raise {
  static_assert(p >= 0, "negative power not allowed in Raise!");
  static constexpr OMEGA_H_INLINE T eval(T x) noexcept {
    return x * Raise<p - 1, T>::eval(x);
  }
};

template <typename T>
struct Raise<0, T> {
  static constexpr OMEGA_H_INLINE T eval(T) noexcept { return 1; }
};

template <Int p, typename T>
constexpr OMEGA_H_INLINE T raise(T x) noexcept {
  return Raise<p, T>::eval(x);
}

template <typename T>
constexpr OMEGA_H_INLINE T square(T x) noexcept {
  return raise<2, T>(x);
}

template <typename T>
OMEGA_H_INLINE T cube(T x) noexcept {
  return raise<3, T>(x);
}

OMEGA_H_INLINE Real sign(Real x) noexcept { return (x < 0.0) ? -1.0 : 1.0; }

OMEGA_H_INLINE Real clamp(Real x, Real low, Real high) noexcept {
  return min2(max2(x, low), high);
}

template <Int dp>
struct Root;

template <>
struct Root<0> {
  static OMEGA_H_INLINE Real eval(Real) noexcept { return 1.0; }
};

template <>
struct Root<1> {
  static OMEGA_H_INLINE Real eval(Real x) noexcept { return x; }
};

template <>
struct Root<2> {
  static OMEGA_H_INLINE Real eval(Real x) noexcept { return std::sqrt(x); }
};

template <>
struct Root<3> {
  static OMEGA_H_INLINE Real eval(Real x) noexcept { return std::cbrt(x); }
};

template <Int p>
OMEGA_H_INLINE Real root(Real x) noexcept {
  return Root<p>::eval(x);
}

/* compile-time-executable Greatest Common Denominator code */
constexpr OMEGA_H_INLINE Int gcd(Int a, Int b) noexcept {
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
  static OMEGA_H_INLINE Real eval(Real x) noexcept {
    return root<dp>(raise<np>(x));
  }
  static_assert(np != dp, "equal case should be specialized");
};

template <Int p>
struct Power<p, p, 1> {
  static OMEGA_H_INLINE Real eval(Real x) noexcept { return x; }
};

template <Int np, Int dp>
OMEGA_H_INLINE Real power(Real x) noexcept {
  return Power<np, dp>::eval(x);
}

OMEGA_H_INLINE Real power(Real x, Int np, Int dp) noexcept {
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
      return -1.0;
    case 2:
      switch (dp) {
        case 1:
          return power<2, 1>(x);
        case 2:
          return power<2, 2>(x);
        case 3:
          return power<2, 3>(x);
      }
      return -1.0;
    case 3:
      switch (dp) {
        case 1:
          return power<3, 1>(x);
        case 2:
          return power<3, 2>(x);
        case 3:
          return power<3, 3>(x);
      }
      return -1.0;
  }
  return -1.0;
}

OMEGA_H_INLINE Real rel_diff_with_floor(
    Real a, Real b, Real floor = EPSILON) noexcept {
  Real am = std::abs(a);
  Real bm = std::abs(b);
  if (am <= floor && bm <= floor) return 0.0;
  return std::abs(b - a) / max2(am, bm);
}

OMEGA_H_INLINE bool are_close(
    Real a, Real b, Real tol = EPSILON, Real floor = EPSILON) noexcept {
  return rel_diff_with_floor(a, b, floor) <= tol;
}

template <typename T>
T divide_no_remainder(T a, T b) {
  OMEGA_H_CHECK(b != 0);
  OMEGA_H_CHECK(a % b == 0);
  return a / b;
}

template <typename T>
struct plus {
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef T result_type;
  OMEGA_H_INLINE T operator()(const T& lhs, const T& rhs) const noexcept {
    return lhs + rhs;
  }
};

template <typename T>
struct logical_and {
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;
  OMEGA_H_INLINE bool operator()(const T& lhs, const T& rhs) const noexcept {
    return lhs && rhs;
  }
};

template <typename T>
struct maximum {
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef T result_type;
  OMEGA_H_INLINE T operator()(const T& lhs, const T& rhs) const noexcept {
    return lhs < rhs ? rhs : lhs;
  }
};

template <typename T>
struct minimum {
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef T result_type;
  OMEGA_H_INLINE T operator()(const T& lhs, const T& rhs) const noexcept {
    return lhs < rhs ? lhs : rhs;
  }
};

template <typename T>
struct identity {
  typedef T argument_type;
  typedef T result_type;
  OMEGA_H_INLINE const T& operator()(const T& x) const noexcept { return x; }
};

template <typename T>
struct multiplies {
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef T result_type;
  OMEGA_H_INLINE T operator()(const T& lhs, const T& rhs) const noexcept {
    return lhs * rhs;
  }
};

// In the algrebra of rotations one often comes across functions that
// take undefined (0/0) values at some points. Close to such points
// these functions must be evaluated using their asymptotic
// expansions; otherwise the computer may produce wildly erroneous
// results or a floating point exception. To avoid unreachable code
// everywhere such functions are used, we introduce here functions to
// the same effect.
//
// Function form: sin(x) / x
// X: 0
// Asymptotics: 1.0 (-x^2/6)
// First radius: (6 * EPS)^(.5)
// Second radius: (120 * EPS)^(.25)
OMEGA_H_INLINE Real sin_x_over_x(Real x) {
  auto const y = std::abs(x);
  auto const e2 = std::sqrt(DBL_EPSILON);
  auto const e4 = std::sqrt(e2);
  if (y > e4) {
    return std::sin(y) / y;
  } else if (y > e2) {
    return 1.0 - y * y / 6.0;
  } else {
    return 1.0;
  }
}

}  // namespace Omega_h

#endif
