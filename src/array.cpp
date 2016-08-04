#include "array.hpp"

#include "algebra.hpp"
#include "functors.hpp"
#include "loop.hpp"

namespace osh {

#ifdef OSH_USE_KOKKOS
template <typename T>
Write<T>::Write(Kokkos::View<T*> view) : view_(view), exists_(true) {}
#endif

template <typename T>
Write<T>::Write(LO size)
    :
#ifdef OSH_USE_KOKKOS
      view_(Kokkos::ViewAllocateWithoutInitializing("omega_h"),
            static_cast<std::size_t>(size))
#else
      ptr_(new T[size], std::default_delete<T[]>()),
      size_(size)
#endif
      ,
      exists_(true) {
}

template <typename T>
static void fill(Write<T> a, T val) {
  auto f = LAMBDA(LO i) { a[i] = val; };
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(LO size, T value) : Write<T>(size) {
  fill(*this, value);
}

template <typename T>
void fill_linear(Write<T> a, T offset, T stride) {
  auto f = LAMBDA(LO i) { a[i] = offset + (stride * static_cast<T>(i)); };
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(LO size, T offset, T stride) : Write<T>(size) {
  fill_linear(*this, offset, stride);
}

template <typename T>
Write<T>::Write(HostWrite<T> host_write) : Write<T>(host_write.write()) {}

template <typename T>
LO Write<T>::size() const {
  CHECK(exists());
#ifdef OSH_USE_KOKKOS
  return static_cast<LO>(view_.size());
#else
  return size_;
#endif
}

/* Several C libraries including ZLib and
   OpenMPI will throw errors when input pointers
   are NULL, even if they point to arrays of size zero. */
template <typename T>
class NonNullPtr {
  static T scratch[1];

 public:
  static T* get(T* p) { return (p == nullptr) ? scratch : p; }
  static T const* get(T const* p) { return (p == nullptr) ? scratch : p; }
};
template <typename T>
T NonNullPtr<T>::scratch[1] = {0};

template <typename T>
T* Write<T>::data() const {
#ifdef OSH_USE_KOKKOS
  return NonNullPtr<T>::get(view_.data());
#else
  return ptr_.get();
#endif
}

#ifdef OSH_USE_KOKKOS
template <typename T>
Kokkos::View<T*> Write<T>::view() const {
  return view_;
}
#endif

template <typename T>
void Write<T>::set(LO i, T value) const {
#ifdef OSH_USE_CUDA
  cudaMemcpy(data() + i, &value, sizeof(T), cudaMemcpyHostToDevice);
#else
  operator[](i) = value;
#endif
}

template <typename T>
T Write<T>::get(LO i) const {
#ifdef OSH_USE_CUDA
  T value;
  cudaMemcpy(&value, data() + i, sizeof(T), cudaMemcpyDeviceToHost);
  return value;
#else
  return operator[](i);
#endif
}

template <typename T>
bool Write<T>::exists() const {
  return exists_;
}

template <typename T>
struct Sum : public SumFunctor<T> {
  using typename SumFunctor<T>::value_type;
  Read<T> a_;
  Sum(Read<T> a) : a_(a) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = update + a_[i];
  }
};

template <typename T>
typename StandinTraits<T>::type sum(Read<T> a) {
  return parallel_reduce(a.size(), Sum<T>(a));
}

template <typename T>
struct Min : public MinFunctor<T> {
  using typename MinFunctor<T>::value_type;
  Read<T> a_;
  Min(Read<T> a) : a_(a) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = min2<value_type>(update, a_[i]);
  }
};

template <typename T>
T min(Read<T> a) {
  auto r = parallel_reduce(a.size(), Min<T>(a));
  return static_cast<T>(r);  // see StandinTraits
}

template <typename T>
struct Max : public MaxFunctor<T> {
  using typename MaxFunctor<T>::value_type;
  Read<T> a_;
  Max(Read<T> a) : a_(a) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = max2<value_type>(update, a_[i]);
  }
};

template <typename T>
T max(Read<T> a) {
  auto r = parallel_reduce(a.size(), Max<T>(a));
  return static_cast<T>(r);  // see StandinTraits
}

template <typename T>
typename StandinTraits<T>::type sum(CommPtr comm, Read<T> a) {
  return comm->allreduce(sum(a), OSH_SUM);
}

template <typename T>
T min(CommPtr comm, Read<T> a) {
  return comm->allreduce(min(a), OSH_MIN);
}

template <typename T>
T max(CommPtr comm, Read<T> a) {
  return comm->allreduce(min(a), OSH_MAX);
}

Reals::Reals() : Read<Real>() {}

Reals::Reals(Write<Real> write) : Read<Real>(write) {}

Reals::Reals(LO size, Real value) : Read<Real>(size, value) {}

Reals::Reals(std::initializer_list<Real> l) : Read<Real>(l) {}

struct AreClose : public AndFunctor {
  Reals a_;
  Reals b_;
  Real tol_;
  Real floor_;
  AreClose(Reals a, Reals b, Real tol, Real floor)
      : a_(a), b_(b), tol_(tol), floor_(floor) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = update && are_close(a_[i], b_[i], tol_, floor_);
  }
};

bool are_close(Reals a, Reals b, Real tol, Real floor) {
  CHECK(a.size() == b.size());
  return static_cast<bool>(
      parallel_reduce(a.size(), AreClose(a, b, tol, floor)));
}

LOs::LOs(Write<LO> write) : Read<LO>(write) {}

LOs::LOs(LO size, LO value) : Read<LO>(size, value) {}

LOs::LOs(LO size, LO offset, LO stride) : Read<LO>(size, offset, stride) {}

LOs::LOs(std::initializer_list<LO> l) : Read<LO>(l) {}

template <typename T>
Read<T>::Read(Write<T> write) : write_(write) {}

template <typename T>
Read<T>::Read(LO size, T value) : write_(size, value) {}

template <typename T>
Read<T>::Read(LO size, T offset, T stride) : write_(size, offset, stride) {}

template <typename T>
Read<T>::Read(std::initializer_list<T> l) : Read<T>(HostWrite<T>(l).write()) {}

template <typename T>
LO Read<T>::size() const {
  return write_.size();
}

template <typename T>
T const* Read<T>::data() const {
  return write_.data();
}

#ifdef OSH_USE_KOKKOS
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

template <typename T>
bool Read<T>::exists() const {
  return write_.exists();
}

template <class T>
struct SameContent : public AndFunctor {
  Read<T> a_;
  Read<T> b_;
  SameContent(Read<T> a, Read<T> b) : a_(a), b_(b) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = update && (a_[i] == b_[i]);
  }
};

template <class T>
bool operator==(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  return parallel_reduce(a.size(), SameContent<T>(a, b));
}

template <class T>
Write<T> deep_copy(Read<T> a) {
  Write<T> b(a.size());
#ifdef OSH_USE_KOKKOS
  Kokkos::deep_copy(b.view(), a.view());
#else
  auto f = LAMBDA(LO i) { b[i] = a[i]; };
  parallel_for(b.size(), f);
#endif
  return b;
}

template <typename T>
HostWrite<T>::HostWrite() = default;

template <typename T>
HostWrite<T>::HostWrite(LO size)
    : write_(size)
#ifdef OSH_USE_KOKKOS
      ,
      mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
}

template <typename T>
HostWrite<T>::HostWrite(LO size, T offset, T stride)
    : HostWrite<T>(Write<T>(size, offset, stride)) {}

template <typename T>
HostWrite<T>::HostWrite(Write<T> write)
    : write_(write)
#ifdef OSH_USE_KOKKOS
      ,
      mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
#ifdef OSH_USE_KOKKOS
  Kokkos::deep_copy(mirror_, write_.view());
#endif
}

template <typename T>
HostWrite<T>::HostWrite(std::initializer_list<T> l)
    :  // an initializer_list should never have over 2 billion items...
      HostWrite<T>(static_cast<LO>(l.size())) {
  LO i = 0;
  for (auto v : l) operator[](i++) = v;
}

template <typename T>
Write<T> HostWrite<T>::write() const {
#ifdef OSH_USE_KOKKOS
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
#ifdef OSH_USE_KOKKOS
  return NonNullPtr<T>::get(mirror_.data());
#else
  return write_.data();
#endif
}

template <typename T>
HostRead<T>::HostRead() = default;

template <typename T>
HostRead<T>::HostRead(Read<T> read)
    : read_(read)
#ifdef OSH_USE_KOKKOS
      ,
      mirror_(Kokkos::create_mirror_view(read.view()))
#endif
{
#ifdef OSH_USE_KOKKOS
  Kokkos::deep_copy(mirror_, read_.view());
#endif
}

template <typename T>
LO HostRead<T>::size() const {
  return read_.size();
}

template <typename T>
T const* HostRead<T>::data() const {
#ifdef OSH_USE_KOKKOS
  return NonNullPtr<T>::get(mirror_.data());
#else
  return read_.data();
#endif
}

template <typename T>
T HostRead<T>::last() const {
  return operator[](size() - 1);
}

template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a) {
  Write<T> b(a.size());
  auto f = LAMBDA(LO i) { b[i] = a[i] * factor; };
  parallel_for(a.size(), f);
  return b;
}

template <typename T>
Read<T> multiply_each(Read<T> a, Read<T> b) {
  CHECK(a.size() % b.size() == 0);
  auto width = a.size() / b.size();
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) {
    for (Int j = 0; j < width; ++j) {
      c[i * width + j] = a[i * width + j] * b[i];
    }
  };
  parallel_for(b.size(), f);
  return c;
}

template <typename T>
Read<T> divide_each(Read<T> a, Read<T> b) {
  CHECK(a.size() % b.size() == 0);
  auto width = a.size() / b.size();
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) {
    for (Int j = 0; j < width; ++j) {
      c[i * width + j] = a[i * width + j] / b[i];
    }
  };
  parallel_for(b.size(), f);
  return c;
}

template <typename T>
Read<T> add_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = a[i] + b[i]; };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> subtract_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = a[i] - b[i]; };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> add_to_each(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = a[i] + b; };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_geq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] >= b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_gt(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] > b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_lt(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] < b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> gt_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] > b[i]); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> geq_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] >= b[i]); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_neq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] != b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_eq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] == b); };
  parallel_for(c.size(), f);
  return c;
}

Read<I8> land_each(Read<I8> a, Read<I8> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] && b[i]); };
  parallel_for(c.size(), f);
  return c;
}

Read<I8> lor_each(Read<I8> a, Read<I8> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] || b[i]); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> get_component(Read<T> a, Int ncomps, Int comp) {
  CHECK(a.size() % ncomps == 0);
  Write<T> b(a.size() / ncomps);
  auto f = LAMBDA(LO i) { b[i] = a[i * ncomps + comp]; };
  parallel_for(b.size(), f);
  return b;
}

#define INST(T)                                                          \
  template class NonNullPtr<T>;                                          \
  template class Write<T>;                                               \
  template class Read<T>;                                                \
  template class HostWrite<T>;                                           \
  template class HostRead<T>;                                            \
  template bool operator==(Read<T> a, Read<T> b);                        \
  template typename StandinTraits<T>::type sum(Read<T> a);               \
  template T min(Read<T> a);                                             \
  template T max(Read<T> a);                                             \
  template typename StandinTraits<T>::type sum(CommPtr comm, Read<T> a); \
  template T min(CommPtr comm, Read<T> a);                               \
  template T max(CommPtr comm, Read<T> a);                               \
  template Write<T> deep_copy(Read<T> a);                                \
  template Read<T> multiply_each_by(T factor, Read<T> x);                \
  template Read<T> multiply_each(Read<T> a, Read<T> b);                  \
  template Read<T> divide_each(Read<T> a, Read<T> b);                    \
  template Read<T> add_each(Read<T> a, Read<T> b);                       \
  template Read<T> subtract_each(Read<T> a, Read<T> b);                  \
  template Read<T> add_to_each(Read<T> a, T b);                          \
  template Read<I8> each_geq_to(Read<T> a, T b);                         \
  template Read<I8> each_gt(Read<T> a, T b);                             \
  template Read<I8> each_lt(Read<T> a, T b);                             \
  template Read<I8> each_neq_to(Read<T> a, T b);                         \
  template Read<I8> each_eq_to(Read<T> a, T b);                          \
  template Read<I8> gt_each(Read<T> a, Read<T> b);                       \
  template Read<T> get_component(Read<T> a, Int ncomps, Int comp);

INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace osh
