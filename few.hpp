template <typename T, Int n>
class Few {
  T array_[n];
  public:
    INLINE T& operator[](Int i) {return array_[i];}
    INLINE T const& operator[](Int i) const {return array_[i];}
    INLINE Few() {}
    Few(std::initializer_list<T> l) {
      Int i = 0;
      for (auto v : l)
        array_[i++] = v;
    }
};
