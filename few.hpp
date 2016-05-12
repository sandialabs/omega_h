template <typename T, U8 n>
class Few {
  T array_[n];
  public:
    INLINE T& operator[](U8 i) {return array_[i];}
};
