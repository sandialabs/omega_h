#ifndef FEW_HPP
#define FEW_HPP

namespace osh {

template <typename T, Int n>
class Few {
  T array_[n];

 public:
  enum { size = n };
  INLINE T& operator[](Int i) { return array_[i]; }
  INLINE T const& operator[](Int i) const { return array_[i]; }
  INLINE Few() {}
  Few(std::initializer_list<T> l) {
    Int i = 0;
    for (auto v : l) array_[i++] = v;
  }
  INLINE void operator=(Few<T, n> const& rhs) volatile {
    for (Int i = 0; i < n; ++i) array_[i] = rhs.array_[i];
  }
  INLINE Few(Few<T, n> const& rhs) {
    for (Int i = 0; i < n; ++i) array_[i] = rhs.array_[i];
  }
  INLINE Few(const volatile Few<T, n>& rhs) {
    for (Int i = 0; i < n; ++i) array_[i] = rhs.array_[i];
  }
};

/* this class avoids dealing with volatile */
template <typename T, Int n>
class HostFew {
  T array_[n];

 public:
  enum { size = n };
  INLINE T& operator[](Int i) { return array_[i]; }
  INLINE T const& operator[](Int i) const { return array_[i]; }
  INLINE HostFew() {}
  HostFew(std::initializer_list<T> l) {
    Int i = 0;
    for (auto v : l) array_[i++] = v;
  }
};

} //end namespace osh

#endif
