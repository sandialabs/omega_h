INLINE Real square(Real x) { return x * x; }
INLINE Real sign(Real x) { return (x < 0) ? -1.0 : 1.0; }

template <U8 n>
class Vector : public Few<Real, n> {
};

/* column-first storage and indexing !!! */
template <U8 m, U8 n>
class Matrix : public Few<Vector<m>, n> {
  public:
    Matrix() {}
    Matrix(std::initializer_list<Real> l) {
      U8 k = 0;
      for (Real v : l) {
        (*this)[k % n][k / n] = v;
        ++k;
      }
    }
};

