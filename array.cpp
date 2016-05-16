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
Write<T>::Write(Int size):
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
  auto f = LAMBDA(Int i) {
    a[i] = val;
  };
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(Int size, T value):
  Write<T>(size)
{
  fill(*this, value);
}

template <typename T>
Int Write<T>::size() const {
#ifdef USE_KOKKOS
  return view_.size();
#else
  return size_;
#endif
}

template <typename T>
T* Write<T>::data() const {
#ifdef USE_KOKKOS
  return view_.data();
#else
  return ptr_.get();
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
  INLINE void operator()(Int i, value_type& update) const
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

Reals::Reals(Int size, Real value):
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
  INLINE void operator()(Int i, value_type& update) const
  {
    update = update && are_close(a_[i], b_[i]);
  }
};

bool are_close(Reals a, Reals b, Real tol, Real floor) {
  CHECK(a.size() == b.size());
  return static_cast<bool>(parallel_reduce(a.size(),
        AreClose(a, b, tol, floor)));
}

LOs::LOs() {
}

LOs::LOs(Read<LO> base):
  Read<LO>(base) {
}

LOs::LOs(Write<LO> write):
  Read<LO>(write) {
}

LOs::LOs(LO size, LO value):
  Read<LO>(size, value) {
}

LOs::LOs(std::initializer_list<LO> l):
  Read<LO>(l) {
}

template <typename T>
Read<T>::Read() {
}

template <typename T>
Read<T>::Read(Write<T> write):
  write_(write) {
}

template <typename T>
Read<T>::Read(Int size, T value):
  write_(size, value) {
}

template <typename T>
Read<T>::Read(std::initializer_list<T> l):
  Read<T>(HostWrite<T>(l).write()) {
}

template <typename T>
Int Read<T>::size() const {
  return write_.size();
}

template <typename T>
T const* Read<T>::data() const {
  return write_.data();
}

#ifdef USE_KOKKOS
template <typename T>
Kokkos::View<const T*> Read<T>::view() const {
  return Kokkos::View<const T*>(write_.view());
}
#endif

template <class T>
struct SameContent : public AndFunctor {
  Read<T> a_;
  Read<T> b_;
  SameContent(Read<T> a, Read<T> b):a_(a),b_(b) {}
  INLINE void operator()(Int i, value_type& update) const
  {
    update = update && (a_[i] == b_[i]);
  }
};

template <class T>
bool operator==(Read<T> a, Read<T> b)
{
  CHECK(a.size() == b.size());
  return parallel_reduce(a.size(), SameContent<T>(a, b));
}

template <typename T>
HostWrite<T>::HostWrite(Int size):
  write_(size)
#ifdef USE_KOKKOS
  ,mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
}

template <typename T>
HostWrite<T>::HostWrite(Write<T> write):
  write_(write)
#ifdef USE_KOKKOS
  ,mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
#ifdef USE_KOKKOS
  Kokkos::deep_copy(mirror_, write_.view());
#endif
}

template <typename T>
HostWrite<T>::HostWrite(std::initializer_list<T> l):
  // an initializer_list should never have over 2 billion items...
  HostWrite<T>(static_cast<Int>(l.size())) {
  Int i = 0;
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
Int HostWrite<T>::size() const {
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
Int HostRead<T>::size() const {
  return read_.size();
}

template <typename T>
Write<T> make_linear(Int n, T offset, T stride) {
  Write<T> a(n);
  auto f = LAMBDA(Int i) {
    a[i] = offset + (stride * static_cast<T>(i));
  };
  parallel_for(n, f);
  return a;
}

#define INST_ARRAY_T(T) \
template class Write<T>; \
template class Read<T>; \
template class HostWrite<T>; \
template class HostRead<T>; \
template bool operator==(Read<T> a, Read<T> b);

INST_ARRAY_T(I8)
INST_ARRAY_T(I16)
INST_ARRAY_T(I32)
INST_ARRAY_T(I64)
INST_ARRAY_T(Real)

template Real sum(Read<Real> a);
template Write<I64> make_linear(Int n, I64 offset, I64 stride);
template Write<I32> make_linear(Int n, I32 offset, I32 stride);
