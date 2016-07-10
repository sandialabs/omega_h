#ifndef TRAITS_HPP
#define TRAITS_HPP

namespace osh {

template <typename T>
struct ArithTraits;

template <>
struct ArithTraits<unsigned char> {
  static INLINE unsigned char max() { return UCHAR_MAX; }
  static INLINE unsigned char min() { return 0; }
};

template <>
struct ArithTraits<signed char> {
  static INLINE signed char max() { return SCHAR_MAX; }
  static INLINE signed char min() { return SCHAR_MIN; }
};

template <>
struct ArithTraits<unsigned int> {
  static INLINE unsigned int max() { return UINT_MAX; }
  static INLINE unsigned int min() { return 0; }
};

template <>
struct ArithTraits<int> {
  static INLINE int max() { return INT_MAX; }
  static INLINE int min() { return INT_MIN; }
};

template <>
struct ArithTraits<unsigned long> {
  static INLINE unsigned long max() { return ULONG_MAX; }
  static INLINE unsigned long min() { return 0; }
};

template <>
struct ArithTraits<signed long> {
  static INLINE signed long max() { return LONG_MAX; }
  static INLINE signed long min() { return LONG_MIN; }
};

template <>
struct ArithTraits<unsigned long long> {
  static INLINE unsigned long long max() { return ULLONG_MAX; }
  static INLINE unsigned long long min() { return 0; }
};

template <>
struct ArithTraits<signed long long> {
  static INLINE signed long long max() { return LLONG_MAX; }
  static INLINE signed long long min() { return LLONG_MIN; }
};

template <>
struct ArithTraits<double> {
  static INLINE double max() { return DBL_MAX; }
  static INLINE double min() { return -DBL_MAX; }
};

} //end namespace osh

#endif
