template <typename T, UInt n>
class Few {
  T array_[n];
  public:
    INLINE T& operator[](UInt i) {return array_[i];}
};
