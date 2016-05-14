struct Int128;

std::ostream& operator<<(std::ostream& o, Int128 const& x);

struct Int128
{
  std::int64_t high;
  std::uint64_t low;
  INLINE Int128() {}
  INLINE Int128(std::int64_t h, std::uint64_t l):
    high(h),low(l) {
  }
  INLINE Int128(std::int64_t value):
    Int128(
      std::int64_t(-1) * (value < 0),
      std::uint64_t(value)) {
  }
  /* volatile... why is this not done by default...
   * returning void instead of reference to *this
   * to silence GCC's warning that the reference
   * is unused. */
  INLINE void operator=(Int128 const& rhs) volatile {
    high = rhs.high;
    low = rhs.low;
  }
  double to_double(int expo) const;
  void print(std::ostream& o) const;
  static INLINE Int128 from_double(double value, int expo) {
    double unit = exp2(double(expo));
    double normalized = value / unit;
    double rounded = round(normalized);
    return Int128(std::int64_t(rounded));
  }
};

/* we moved the actual operators out here to take
 * their arguments by value, the CUDA compiler
 * wouldn't match them otherwise when operating on
 * two volatile Int128 variables.
 */

INLINE Int128 operator+(Int128 lhs, Int128 rhs) {
  Int128 sum;
  sum.high = lhs.high + rhs.high;
  sum.low = lhs.low + rhs.low;
  // check for overflow of low 64 bits, add carry to high
  sum.high += (sum.low < lhs.low);
  return sum;
}

INLINE Int128 operator-(Int128 lhs, Int128 rhs) {
  Int128 difference;
  difference.high = lhs.high - rhs.high;
  difference.low = lhs.low - rhs.low;
  // check for underflow of low 64 bits, subtract carry from high
  difference.high -= (difference.low > lhs.low);
  return difference;
}

INLINE Int128 operator-(Int128 x) {
  return Int128(0) - x;
}

INLINE Int128 operator>>(Int128 x, int expo) {
  Int128 shifted;
  shifted.low = (x.low >> expo) |
    (std::uint64_t(x.high) << (64 - expo));
  shifted.high = x.high >> expo;
  return shifted;
}

INLINE bool operator==(Int128 lhs, Int128 rhs) {
  return lhs.high == rhs.high && lhs.low == rhs.low;
}

INLINE bool operator<(Int128 lhs, Int128 rhs) {
  if (lhs.high != rhs.high)
    return lhs.high < rhs.high;
  return lhs.low < rhs.low;
}
