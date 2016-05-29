template <typename T>
class HostWrite;

template <typename T>
class Write {
#ifdef USE_KOKKOS
  Kokkos::View<T*> view_;
#else
  std::shared_ptr<T> ptr_;
  LO size_;
#endif
public:
  Write();
  Write(LO size);
  Write(LO size, T value);
  Write(LO size, T offset, T stride);
  Write(HostWrite<T> host_write);
  LO size() const;
  INLINE T& operator[](LO i) const {
#ifdef CHECK_BOUNDS
    if (i < 0)
      std::cerr << "i = " << i << '\n';
    CHECK(0 <= i);
    if (i >= size())
      std::cerr << "i = " << i << '\n';
    CHECK(i < size());
#endif
#ifdef USE_KOKKOS
    return view_(i);
#else
    return ptr_.get()[i];
#endif
  }
  T* data() const;
#ifdef USE_KOKKOS
  Kokkos::View<T*> view() const;
#endif
  void set(LO i, T value) const;
  T get(LO i) const;
  bool exists() const;
};

template <typename T>
class Read {
  Write<T> write_;
public:
  Read();
  Read(Write<T> write);
  Read(LO size, T value);
  Read(LO size, T offset, T stride);
  Read(std::initializer_list<T> l);
  LO size() const;
  INLINE T const& operator[](LO i) const {
    return write_[i];
  }
  T const* data() const;
#ifdef USE_KOKKOS
  Kokkos::View<const T*> view() const;
#endif
  T get(LO i) const;
  T last() const;
  bool exists() const;
};

template <class T>
bool operator==(Read<T> a, Read<T> b);

template <class T>
Write<T> deep_copy(Read<T> a);

template <typename T>
T sum(Read<T> a);
template <typename T>
T max(Read<T> a);

class Reals : public Read<Real> {
public:
  Reals();
  Reals(Read<Real> base);
  Reals(Write<Real> write);
  Reals(LO size, Real value);
  Reals(std::initializer_list<Real> l);
};

bool are_close(Reals a, Reals b, Real tol = EPSILON, Real floor = EPSILON);

class LOs : public Read<LO> {
public:
  LOs();
  LOs(Read<LO> base);
  LOs(Write<LO> write);
  LOs(LO size, LO value);
  LOs(LO size, LO offset, LO stride);
  LOs(std::initializer_list<LO> l);
};

template <typename T>
class HostWrite {
  Write<T> write_;
#ifdef USE_KOKKOS
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
#ifdef USE_KOKKOS
#ifdef CHECK_BOUNDS
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
class HostRead {
  Read<T> read_;
#ifdef USE_KOKKOS
  typename Kokkos::View<const T*>::HostMirror mirror_;
#endif
public:
  HostRead();
  HostRead(Read<T> read);
  LO size() const;
  inline T const& operator[](LO i) const {
#ifdef USE_KOKKOS
#ifdef CHECK_BOUNDS
    CHECK(0 <= i);
    CHECK(i < size());
#endif
    return mirror_(i);
#else
    return read_[i];
#endif
  }
  T const* data() const;
  T last() const;
};

template <typename T>
std::ostream& operator<<(std::ostream& o, Read<T> a);

template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a);
template <typename T>
Read<T> multiply_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_each(Read<T> a, Read<T> b);
