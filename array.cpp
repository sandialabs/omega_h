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
Write<T>::Write(LO size):
#ifdef USE_KOKKOS
  view_("omega_h", static_cast<std::size_t>(size))
#else
  ptr_(new T[size], std::default_delete<T[]>()),
  size_(size)
#endif
{
}

template <typename T>
static void fill(Write<T> a, T val) {
  auto f = LAMBDA(LO i) {
    a[i] = val;
  };
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(LO size, T value):
  Write<T>(size)
{
  fill(*this, value);
}

template <typename T>
void fill_linear(Write<T> a, T offset, T stride) {
  auto f = LAMBDA(LO i) {
    a[i] = offset + (stride * static_cast<T>(i));
  };
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(LO size, T offset, T stride):
  Write<T>(size)
{
  fill_linear(*this, offset, stride);
}

template <typename T>
Write<T>::Write(HostWrite<T> host_write):
  Write<T>(host_write.write())
{
}

template <typename T>
LO Write<T>::size() const {
#ifdef USE_KOKKOS
  return static_cast<LO>(view_.size());
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
void Write<T>::set(LO i, T value) const {
#ifdef USE_CUDA
  cudaMemcpy(data() + i, &value, sizeof(T), cudaMemcpyHostToDevice);
#else
  operator[](i) = value;
#endif
}

template <typename T>
T Write<T>::get(LO i) const {
#ifdef USE_CUDA
  T value;
  cudaMemcpy(&value, data() + i, sizeof(T), cudaMemcpyDeviceToHost);
  return value;
#else
  return operator[](i);
#endif
}

template <typename T>
struct Sum : public SumFunctor<T> {
  typedef typename SumFunctor<T>::value_type value_type;
  Read<T> a_;
  Sum(Read<T> a):a_(a) {}
  INLINE void operator()(LO i, value_type& update) const
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

Reals::Reals(Read<Real> base):
  Read<Real>(base)
{}

Reals::Reals(Write<Real> write):
  Read<Real>(write)
{}

Reals::Reals(LO size, Real value):
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
  INLINE void operator()(LO i, value_type& update) const
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

LOs::LOs(LO size, LO offset, LO stride):
  Read<LO>(size, offset, stride) {
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
Read<T>::Read(LO size, T value):
  write_(size, value) {
}

template <typename T>
Read<T>::Read(LO size, T offset, T stride):
  write_(size, offset, stride) {
}

template <typename T>
Read<T>::Read(std::initializer_list<T> l):
  Read<T>(HostWrite<T>(l).write()) {
}

template <typename T>
LO Read<T>::size() const {
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

template <typename T>
T Read<T>::get(LO i) const {
  return write_.get(i);
}

template <typename T>
T Read<T>::last() const {
  return get(size() - 1);
}

template <class T>
struct SameContent : public AndFunctor {
  Read<T> a_;
  Read<T> b_;
  SameContent(Read<T> a, Read<T> b):a_(a),b_(b) {}
  INLINE void operator()(LO i, value_type& update) const
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

template <class T>
Write<T> deep_copy(Read<T> a) {
  Write<T> b(a.size());
  auto f = LAMBDA(LO i) {
    b[i] = a[i];
  };
  parallel_for(b.size(), f);
  return b;
}

template <typename T>
HostWrite<T>::HostWrite(LO size):
  write_(size)
#ifdef USE_KOKKOS
  ,mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
}

template <typename T>
HostWrite<T>::HostWrite(LO size, T offset, T stride):
  HostWrite<T>(Write<T>(size, offset, stride))
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
  HostWrite<T>(static_cast<LO>(l.size())) {
  LO i = 0;
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
LO HostWrite<T>::size() const {
  return write_.size();
}

template <typename T>
T* HostWrite<T>::data() const {
#ifdef USE_KOKKOS
  return mirror_.data();
#else
  return write_.data();
#endif
}

template <typename T>
HostRead<T>::HostRead() {
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
LO HostRead<T>::size() const {
  return read_.size();
}

template <typename T>
T const* HostRead<T>::data() const {
#ifdef USE_KOKKOS
  return mirror_.data();
#else
  return read_.data();
#endif
}

template <typename T>
std::ostream& operator<<(std::ostream& o, Read<T> a) {
  HostRead<T> ha = a;
  if (ha.size() <= 20) {
    for (LO i = 0; i < ha.size(); ++i)
      o << ' ' << ha[i];
    o << '\n';
  } else {
    for (LO i = 0; i < ha.size(); ++i)
      o << ha[i] << '\n';
  }
  return o;
}

template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a) {
  Write<T> b(a.size());
  auto f = LAMBDA(LO i) {
    b[i] = a[i] * factor;
  };
  parallel_for(a.size(), f);
  return b;
}

template <typename T>
Read<T> add_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) {
    c[i] = a[i] + b[i];
  };
  parallel_for(c.size(), f);
  return c;
}

#define INST_ARRAY_T(T) \
template class Write<T>; \
template class Read<T>; \
template class HostWrite<T>; \
template class HostRead<T>; \
template bool operator==(Read<T> a, Read<T> b); \
template Write<T> deep_copy(Read<T> a); \
template std::ostream& operator<<(std::ostream& o, Read<T> a); \
template Read<T> multiply_each_by(T factor, Read<T> x); \
template Read<T> add_each(Read<T> a, Read<T> b);

INST_ARRAY_T(I8)
INST_ARRAY_T(I16)
INST_ARRAY_T(I32)
INST_ARRAY_T(I64)
INST_ARRAY_T(Real)

template Real sum(Read<Real> a);
