template <typename T, UInt n>
class Few {
  T array_[n];
  public:
    INLINE T& operator[](UInt i) {return array_[i];}
    INLINE Few() {}
    Few(std::initializer_list<T> l) {
      UInt i = 0;
      for (auto v : l)
        array_[i++] = v;
    }
};
