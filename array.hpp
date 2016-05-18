template <typename T>
class HostWrite;

template <typename T>
class Write {
#ifdef USE_KOKKOS
  Kokkos::View<T*> view_;
#else
  std::shared_ptr<T> ptr_;
  Int size_;
#endif
public:
  Write();
  Write(Int size);
  Write(Int size, T value);
  Write(Int size, T offset, T stride);
  Write(HostWrite<T> host_write);
  Int size() const;
  INLINE T& operator[](Int i) const {
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
  void set(Int i, T value) const;
  T get(Int i) const;
};

template <typename T>
class Read {
  Write<T> write_;
public:
  Read();
  Read(Write<T> write);
  Read(Int size, T value);
  Read(Int size, T offset, T stride);
  Read(std::initializer_list<T> l);
  Int size() const;
  INLINE T const& operator[](Int i) const {
    return write_[i];
  }
  T const* data() const;
#ifdef USE_KOKKOS
  Kokkos::View<const T*> view() const;
#endif
  T get(Int i) const;
};

template <class T>
bool operator==(Read<T> a, Read<T> b);

template <typename T>
T sum(Read<T> a);

class Reals : public Read<Real> {
public:
  Reals();
  Reals(Read<Real> base);
  Reals(Write<Real> write);
  Reals(Int size, Real value);
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
  HostWrite(Int size);
  HostWrite(Int size, T offset, T stride);
  HostWrite(Write<T> write);
  HostWrite(std::initializer_list<T> l);
  Write<T> write() const;
  Int size() const;
  inline T& operator[](Int i) const {
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
};

template <typename T>
class HostRead {
  Read<T> read_;
#ifdef USE_KOKKOS
  typename Kokkos::View<const T*>::HostMirror mirror_;
#endif
public:
  HostRead(Read<T> read);
  Int size() const;
  inline T const& operator[](Int i) const {
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
};

template <typename T>
std::ostream& operator<<(std::ostream& o, Read<T> a);
