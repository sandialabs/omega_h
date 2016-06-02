template <typename T, Int n>
class Few {
  T array_[n];
  public:
    OSH_INLINE T& operator[](Int i) {return array_[i];}
    OSH_INLINE T const& operator[](Int i) const {return array_[i];}
    OSH_INLINE Few() {}
    Few(std::initializer_list<T> l) {
      Int i = 0;
      for (auto v : l)
        array_[i++] = v;
    }
    OSH_INLINE void operator=(Few<T,n> const& rhs) volatile {
      for (Int i = 0; i < n; ++i)
        array_[i] = rhs.array_[i];
    }
    OSH_INLINE Few(Few<T,n> const& rhs) {
      for (Int i = 0; i < n; ++i)
        array_[i] = rhs.array_[i];
    }
    OSH_INLINE Few(const volatile Few<T,n>& rhs) {
      for (Int i = 0; i < n; ++i)
        array_[i] = rhs.array_[i];
    }
};
