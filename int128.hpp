class Int128
{
  std::uint64_t low;
  std::int64_t high;
public:
  INLINE Int128() {}
  INLINE Int128(std::int64_t value):
    low(std::uint64_t(value)),
    high(std::int64_t(-1) * (value < 0)) {
  }
  INLINE Int128(double value, int expo):
    Int128(std::int64_t(round(value / exp2(double(expo))))) {
  }
  INLINE Int128 operator+(const Int128 & rhs) {
    Int128 sum;
    sum.high = high + rhs.high;
    sum.low = low + rhs.low;
    // check for overflow of low 64 bits, add carry to high
    sum.high += (sum.low < low);
    return sum;
  }
  INLINE Int128 operator-(const Int128 & rhs) {
    Int128 difference;
    difference.high = high - rhs.high;
    difference.low = low - rhs.low;
    // check for underflow of low 64 bits, subtract carry to high
    difference.high -= (difference.low > low);
    return difference;
  }
  INLINE Int128 operator>>(int expo) {
    Int128 shifted;
    shifted.low = (low >> expo) |
      (std::uint64_t(high) << (64 - expo));
    shifted.high = high >> expo;
    return shifted;
  }
  void print(std::ostream& o);
};

std::ostream& operator<<(std::ostream& o, Int128 x);
