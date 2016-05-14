template <typename T>
Write<T>::Write():
#ifdef USE_KOKKOS
  view_()
#else
  ptr_(),
  size_(0)
#endif
{
}

template <typename T>
Write<T>::Write(UInt size):
#ifdef USE_KOKKOS
  view_("omega_h", size)
#else
  ptr_(new T[size], std::default_delete<T[]>()),
  size_(size)
#endif
{
}

template <typename T>
static void fill(Write<T> a, T val) {
  auto f = LAMBDA(UInt i) {
    a[i] = val;
  };
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(UInt size, T value):
  Write<T>(size)
{
  fill(*this, value);
}

template <typename T>
UInt Write<T>::size() const {
#ifdef USE_KOKKOS
  return view_.size();
#else
  return size_;
#endif
}

#ifdef USE_KOKKOS
template <typename T>
Kokkos::View<T*> Write<T>::view() const {
  return view_;
}
#endif

template <typename T>
struct Sum : public SumFunctor<T> {
  typedef typename SumFunctor<T>::value_type value_type;
  Read<T> a_;
  Sum(Read<T> a):a_(a) {}
  INLINE void operator()(UInt i, value_type& update) const
  {
    update = update + a_[i];
  }
};

template <typename T>
T sum(Read<T> a) {
  return parallel_reduce(a.size(), Sum<T>(a));
}

Reals::Reals():
  Read<Real>()
{}

Reals::Reals(Write<Real> write):
  Read<Real>(write)
{}

Reals::Reals(UInt size, Real value):
  Read<Real>(size, value)
{}

Reals::Reals(std::initializer_list<Real> l):
  Read<Real>(l)
{}

struct AreClose : public AndFunctor {
  Reals a_;
  Reals b_;
  Real tol_;
  Real floor_;
  AreClose(Reals a, Reals b, Real tol, Real floor):
    a_(a),b_(b),tol_(tol),floor_(floor)
  {}
  INLINE void operator()(UInt i, value_type& update) const
  {
    update = update && are_close(a_[i], b_[i]);
  }
};

bool are_close(Reals a, Reals b, Real tol, Real floor) {
  CHECK(a.size() == b.size());
  return static_cast<bool>(parallel_reduce(a.size(),
        AreClose(a, b, tol, floor)));
}

template <typename T>
Read<T>::Read() {
}

template <typename T>
Read<T>::Read(Write<T> write):
  write_(write) {
}

template <typename T>
Read<T>::Read(UInt size, T value):
  write_(size, value) {
}

template <typename T>
Read<T>::Read(std::initializer_list<T> l):
  Read<T>(HostWrite<T>(l).write()) {
}

template <typename T>
UInt Read<T>::size() const {
  return write_.size();
}

#ifdef USE_KOKKOS
template <typename T>
Kokkos::View<const T*> Read<T>::view() const {
  return Kokkos::View<const T*>(write_.view());
}
#endif

template <typename T>
HostWrite<T>::HostWrite(UInt size):
  write_(size)
#ifdef USE_KOKKOS
  ,mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
}

template <typename T>
HostWrite<T>::HostWrite(std::initializer_list<T> l):
  HostWrite<T>(l.size()) {
  UInt i = 0;
  for (auto v : l)
    operator[](i++) = v;
}

template <typename T>
Write<T> HostWrite<T>::write() const {
#ifdef USE_KOKKOS
  Kokkos::deep_copy(write_.view(), mirror_);
#endif
  return write_;
}

template <typename T>
UInt HostWrite<T>::size() const {
  return write_.size();
}

template <typename T>
HostRead<T>::HostRead(Read<T> read):
  read_(read)
#ifdef USE_KOKKOS
  ,mirror_(Kokkos::create_mirror_view(read.view()))
#endif
{
#ifdef USE_KOKKOS
  Kokkos::deep_copy(mirror_, read_.view());
#endif
}

template <typename T>
UInt HostRead<T>::size() const {
  return read_.size();
}

#define INST_ARRAY_T(T) \
template class Write<T>; \
template class Read<T>; \
template class HostWrite<T>; \
template class HostRead<T>;

INST_ARRAY_T(U8)
INST_ARRAY_T(U16)
INST_ARRAY_T(U32)
INST_ARRAY_T(U64)
INST_ARRAY_T(Real)

template Real sum(Read<Real> a);
