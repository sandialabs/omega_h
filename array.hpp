template <typename T>
class HostWrite;

template <class T>
bool operator==(Read<T> a, Read<T> b);

template <class T>
Write<T> deep_copy(Read<T> a);

template <typename T>
typename StandinTraits<T>::type sum(Read<T> a);
template <typename T>
T min(Read<T> a);
template <typename T>
T max(Read<T> a);

bool are_close(Reals a, Reals b, Real tol = EPSILON, Real floor = EPSILON);

template <typename T>
class HostWrite {
  Write<T> write_;
#ifdef OSH_USE_KOKKOS
  typename Kokkos::View<T*>::HostMirror mirror_;
#endif
public:
  HostWrite(LO size);
  HostWrite(LO size, T offset, T stride);
  HostWrite(Write<T> write);
  HostWrite(std::initializer_list<T> l);
  Write<T> write() const;
  LO size() const;
  inline T& operator[](LO i) const {
#ifdef OSH_USE_KOKKOS
#ifdef OSH_CHECK_BOUNDS
    CHECK(0 <= i);
    CHECK(i < size());
#endif
    return mirror_(i);
#else
    return write_[i];
#endif
  }
  T* data() const;
};

template <typename T>
std::ostream& operator<<(std::ostream& o, Read<T> a);

template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a);
template <typename T>
Read<T> multiply_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> divide_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> subtract_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_to_each(Read<T> a, T b);
template <typename T>
Read<I8> each_geq_to(Read<T> a, T b);
template <typename T>
Read<I8> each_gt(Read<T> a, T b);
template <typename T>
Read<I8> each_lt(Read<T> a, T b);
template <typename T>
Read<I8> gt_each(Read<T> a, Read<T> b);
template <typename T>
Read<I8> each_neq_to(Read<T> a, T b);
template <typename T>
Read<I8> each_eq_to(Read<T> a, T b);
Read<I8> land_each(Read<I8> a, Read<I8> b);
Read<I8> lor_each(Read<I8> a, Read<I8> b);
