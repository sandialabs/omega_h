#ifndef OMEGA_H_INT128_HPP
#define OMEGA_H_INT128_HPP

#include <cstdint>
#include <iosfwd>

#include <Omega_h_macros.h>

namespace Omega_h {

struct Int128 {
  std::int64_t high;
  std::uint64_t low;
  OMEGA_H_INLINE Int128();
  OMEGA_H_INLINE Int128(std::int64_t h, std::uint64_t l);
  OMEGA_H_INLINE Int128(std::int64_t value);
  OMEGA_H_INLINE void operator=(Int128 const& rhs) volatile;
  OMEGA_H_INLINE Int128(Int128 const& rhs);
  OMEGA_H_INLINE Int128(const volatile Int128& rhs);
  double to_double(double unit) const;
  void print(std::ostream& o) const;
  static OMEGA_H_INLINE Int128 from_double(double value, double unit);
};

/*
   We code our own int128 because we target GPUs and similar
   systems where such a type is not guaranteed to exist.
   TODO: wrap built-in types when they are available
*/

OMEGA_H_INLINE Int128::Int128() {}

OMEGA_H_INLINE Int128::Int128(std::int64_t h, std::uint64_t l)
    : high(h), low(l) {}

OMEGA_H_INLINE Int128::Int128(std::int64_t value)
    : Int128(std::int64_t(-1) * (value < 0), std::uint64_t(value)) {}

/* volatile... why is this not done by default...
 * returning void instead of reference to *this
 * to silence GCC's warning that the reference
 * is unused.
 */
OMEGA_H_INLINE void Int128::operator=(Int128 const& rhs) volatile {
  high = rhs.high;
  low = rhs.low;
}

/* which implies we have to declare a regular copy
   constructor */
OMEGA_H_INLINE Int128::Int128(Int128 const& rhs)
    : high(rhs.high), low(rhs.low) {}

/* and a volatile rhs one ? */
OMEGA_H_INLINE Int128::Int128(const volatile Int128& rhs)
    : high(rhs.high), low(rhs.low) {}

OMEGA_H_INLINE Int128 Int128::from_double(double value, double unit) {
  double normalized = value / unit;
  return Int128(std::int64_t(normalized));
}

/* we moved the actual operators out here to take
 * their arguments by value, the CUDA compiler
 * wouldn't match them otherwise when operating on
 * two volatile Int128 variables.
 */

OMEGA_H_INLINE Int128 operator+(Int128 lhs, Int128 rhs) {
  Int128 sum;
  sum.high = lhs.high + rhs.high;
  sum.low = lhs.low + rhs.low;
  // check for overflow of low 64 bits, add carry to high
  sum.high += (sum.low < lhs.low);
  return sum;
}

OMEGA_H_INLINE Int128 operator-(Int128 lhs, Int128 rhs) {
  Int128 difference;
  difference.high = lhs.high - rhs.high;
  difference.low = lhs.low - rhs.low;
  // check for underflow of low 64 bits, subtract carry from high
  difference.high -= (difference.low > lhs.low);
  return difference;
}

OMEGA_H_INLINE Int128 operator-(Int128 x) { return Int128(0) - x; }

OMEGA_H_INLINE Int128 operator>>(Int128 x, int expo) {
  Int128 shifted;
  shifted.low = (x.low >> expo) | (std::uint64_t(x.high) << (64 - expo));
  shifted.high = x.high >> expo;
  return shifted;
}

OMEGA_H_INLINE bool operator==(Int128 lhs, Int128 rhs) {
  return lhs.high == rhs.high && lhs.low == rhs.low;
}

OMEGA_H_INLINE bool operator<(Int128 lhs, Int128 rhs) {
  if (lhs.high != rhs.high) return lhs.high < rhs.high;
  return lhs.low < rhs.low;
}

}  // namespace Omega_h

#endif
