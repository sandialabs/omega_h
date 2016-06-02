template <typename T>
OSH_INLINE T max2(T a, T b) { return (b > a) ? (b) : (a); }
template <typename T>
OSH_INLINE T min2(T a, T b) { return (b < a) ? (b) : (a); }
template <typename T>
OSH_INLINE void swap2(T& a, T& b) {
  T c = a;
  a = b;
  b = c;
}
