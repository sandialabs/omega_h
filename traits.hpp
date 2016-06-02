template <typename T>
struct ArithTraits;

template <>
struct ArithTraits<unsigned char> {
  static OSH_INLINE unsigned char max() { return UCHAR_MAX; }
  static OSH_INLINE unsigned char min() { return 0; }
};

template <>
struct ArithTraits<signed char> {
  static OSH_INLINE signed char max() { return SCHAR_MAX; }
  static OSH_INLINE signed char min() { return SCHAR_MIN; }
};

template <>
struct ArithTraits<unsigned int> {
  static OSH_INLINE unsigned int max() { return UINT_MAX; }
  static OSH_INLINE unsigned int min() { return 0; }
};

template <>
struct ArithTraits<int> {
  static OSH_INLINE int max() { return INT_MAX; }
  static OSH_INLINE int min() { return INT_MIN; }
};

template <>
struct ArithTraits<unsigned long> {
  static OSH_INLINE unsigned long max() { return ULONG_MAX; }
  static OSH_INLINE unsigned long min() { return 0; }
};

template <>
struct ArithTraits<signed long> {
  static OSH_INLINE signed long max() { return LONG_MAX; }
  static OSH_INLINE signed long min() { return LONG_MIN; }
};

template <>
struct ArithTraits<unsigned long long> {
  static OSH_INLINE unsigned long long max() { return ULLONG_MAX; }
  static OSH_INLINE unsigned long long min() { return 0; }
};

template <>
struct ArithTraits<signed long long> {
  static OSH_INLINE signed long long max() { return LLONG_MAX; }
  static OSH_INLINE signed long long min() { return LLONG_MIN; }
};

template <>
struct ArithTraits<double> {
  static OSH_INLINE double max() { return DBL_MAX; }
  static OSH_INLINE double min() { return -DBL_MAX; }
};
