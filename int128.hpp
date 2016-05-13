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
  INLINE Int128 operator+(Int128 const& rhs) const {
    Int128 sum;
    sum.high = high + rhs.high;
    sum.low = low + rhs.low;
    // check for overflow of low 64 bits, add carry to high
    sum.high += (sum.low < low);
    return sum;
  }
  INLINE Int128 operator-(Int128 const& rhs) const {
    Int128 difference;
    difference.high = high - rhs.high;
    difference.low = low - rhs.low;
    // check for underflow of low 64 bits, subtract carry from high
    difference.high -= (difference.low > low);
    return difference;
  }
  INLINE Int128 operator-() const {
    return Int128(0) - *this;
  }
  INLINE Int128 operator>>(int expo) const {
    Int128 shifted;
    shifted.low = (low >> expo) |
      (std::uint64_t(high) << (64 - expo));
    shifted.high = high >> expo;
    return shifted;
  }
  INLINE bool operator==(Int128 const& rhs) {
    return high == rhs.high && low == rhs.low;
  }
  INLINE bool operator<(Int128 const& rhs) {
    if (high != rhs.high)
      return high < rhs.high;
    return low < rhs.low;
  }
  INLINE double as_double(int expo) {
    Int128 tmp = *this;
    if (tmp < Int128(0))
      tmp = -tmp;
    while (tmp.high) {
      tmp = tmp >> 1;
      ++expo;
    }
    double x = tmp.low;
    if (*this < Int128(0))
      x = -x;
    x *= exp2(double(expo));
    return x;
  }
  static INLINE Int128 max() {
    Int128 x;
    x.high = INT64_MAX;
    x.low = UINT64_MAX;
    return x;
  }
  static INLINE Int128 min() {
    Int128 x;
    x.high = INT64_MIN;
    x.low = 0;
    return x;
  }
  void print(std::ostream& o);
};

std::ostream& operator<<(std::ostream& o, Int128 x);
