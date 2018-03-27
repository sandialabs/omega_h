#include "Omega_h_array.hpp"

#include <cstring>
#include <sstream>

#include "Omega_h_control.hpp"
#include "Omega_h_functors.hpp"
#include "Omega_h_loop.hpp"

namespace Omega_h {

/* Several C libraries including ZLib and
   OpenMPI will throw errors when input pointers
   are NULL, even if they point to arrays of size zero. */
template <typename T>
class NonNullPtr {
  static T scratch[1];

 public:
  static T* get(T* p) { return (p == nullptr) ? scratch : p; }
};
template <typename T>
T NonNullPtr<T>::scratch[1] = {0};

template <typename T>
T* nonnull(T* p) {
  return NonNullPtr<T>::get(p);
}

static std::size_t current_array_bytes = 0;

std::size_t get_current_bytes() { return current_array_bytes; }

static std::size_t max_array_bytes = 0;

std::size_t get_max_bytes() { return max_array_bytes; }

template <typename T>
void Write<T>::log_allocation() const {
  if (!should_log_memory) return;
  current_array_bytes += bytes();
  if (current_array_bytes > max_array_bytes) {
    max_array_bytes = current_array_bytes;
    delete[] max_memory_stacktrace;
    max_memory_stacktrace = nullptr;
    std::stringstream ss;
    print_stacktrace(ss, 64);
    auto s = ss.str();
    max_memory_stacktrace = new char[s.length() + 1];
    strcpy(max_memory_stacktrace, s.c_str());
  }
}

#ifdef OMEGA_H_USE_KOKKOSCORE
template <typename T>
Write<T>::Write(Kokkos::View<T*> view_in) : view_(view_in) {
  log_allocation();
}
#endif

template <typename T>
Write<T>::Write(LO size_in, std::string const& name) {
// begin_code("Write(size,name)");
#ifdef OMEGA_H_USE_KOKKOSCORE
  view_ = decltype(view_)(Kokkos::ViewAllocateWithoutInitializing(name),
      static_cast<std::size_t>(size_in));
#else
  (void)name;
  ptr_ = decltype(ptr_)(new T[size_in], std::default_delete<T[]>());
  size_ = size_in;
#endif
  log_allocation();
  // end_code();
}

template <typename T>
void Write<T>::check_release() const {
  if (should_log_memory && use_count() == 1) {
    current_array_bytes -= bytes();
  }
}

template <typename T>
Write<T>& Write<T>::operator=(Write<T> const& other) {
  check_release();
#ifdef OMEGA_H_USE_KOKKOSCORE
  view_ = other.view_;
#else
  ptr_ = other.ptr_;
  size_ = other.size_;
#endif
  return *this;
}

template <typename T>
static void fill(Write<T> a, T val) {
  auto f = OMEGA_H_LAMBDA(LO i) { a[i] = val; };
  parallel_for(a.size(), f, "Write(size,value)");
}

template <typename T>
Write<T>::Write(LO size_in, T value, std::string const& name)
    : Write<T>(size_in, name) {
  fill(*this, value);
}

template <typename T>
void fill_linear(Write<T> a, T offset, T stride) {
  auto f = OMEGA_H_LAMBDA(LO i) {
    a[i] = offset + (stride * static_cast<T>(i));
  };
  parallel_for(a.size(), f, "Write(size,offset,stride)");
}

template <typename T>
Write<T>::Write(LO size_in, T offset, T stride, std::string const& name)
    : Write<T>(size_in, name) {
  fill_linear(*this, offset, stride);
}

template <typename T>
Write<T>::Write(HostWrite<T> host_write) : Write<T>(host_write.write()) {}

template <typename T>
std::size_t Write<T>::bytes() const {
  return static_cast<std::size_t>(size()) * sizeof(T);
}

template <typename T>
void Write<T>::set(LO i, T value) const {
#ifdef OMEGA_H_USE_CUDA
  cudaMemcpy(data() + i, &value, sizeof(T), cudaMemcpyHostToDevice);
#else
  operator[](i) = value;
#endif
}

template <typename T>
T Write<T>::get(LO i) const {
#ifdef OMEGA_H_USE_CUDA
  T value;
  cudaMemcpy(&value, data() + i, sizeof(T), cudaMemcpyDeviceToHost);
  return value;
#else
  return operator[](i);
#endif
}

Bytes::Bytes(Write<Byte> write) : Read<Byte>(write) {}

Bytes::Bytes(LO size_in, Byte value, std::string const& name)
    : Read<Byte>(size_in, value, name) {}

Bytes::Bytes(std::initializer_list<Byte> l, std::string const& name)
    : Read<Byte>(l, name) {}

LOs::LOs(Write<LO> write) : Read<LO>(write) {}

LOs::LOs(LO size_in, LO value, std::string const& name)
    : Read<LO>(size_in, value, name) {}

LOs::LOs(LO size_in, LO offset, LO stride, std::string const& name)
    : Read<LO>(size_in, offset, stride, name) {}

LOs::LOs(std::initializer_list<LO> l, std::string const& name)
    : Read<LO>(l, name) {}

GOs::GOs(Write<GO> write) : Read<GO>(write) {}

GOs::GOs(LO size_in, GO value, std::string const& name)
    : Read<GO>(size_in, value, name) {}

GOs::GOs(LO size_in, GO offset, GO stride, std::string const& name)
    : Read<GO>(size_in, offset, stride, name) {}

GOs::GOs(std::initializer_list<GO> l, std::string const& name)
    : Read<GO>(l, name) {}

Reals::Reals(Write<Real> write) : Read<Real>(write) {}

Reals::Reals(LO size_in, Real value, std::string const& name)
    : Read<Real>(size_in, value, name) {}

Reals::Reals(std::initializer_list<Real> l, std::string const& name)
    : Read<Real>(l, name) {}

template <typename T>
Read<T>::Read(Write<T> write) : write_(write) {}

template <typename T>
Read<T>::Read(LO size_in, T value, std::string const& name)
    : Read<T>(Write<T>(size_in, value, name)) {}

template <typename T>
Read<T>::Read(LO size_in, T offset, T stride, std::string const& name)
    : Read<T>(Write<T>(size_in, offset, stride, name)) {}

template <typename T>
Read<T>::Read(std::initializer_list<T> l, std::string const& name)
    : Read<T>(HostWrite<T>(l, name).write()) {}

#ifdef OMEGA_H_USE_KOKKOSCORE
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
T Read<T>::first() const {
  return get(0);
}

template <typename T>
T Read<T>::last() const {
  return get(size() - 1);
}

template <typename T>
HostWrite<T>::HostWrite() = default;

#ifdef OMEGA_H_USE_KOKKOSCORE
template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_uninit_mirror(
    const Kokkos::View<T, P...>& src) {
  typedef Kokkos::View<T, P...> src_type;
  typedef typename src_type::HostMirror dst_type;
  return dst_type(Kokkos::ViewAllocateWithoutInitializing(
                      std::string("host_") + src.label()),
      src.extent(0), src.extent(1), src.extent(2),
      src.extent(3), src.extent(4), src.extent(5),
      src.extent(6), src.extent(7));
}

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_uninit_mirror_view(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<(
        std::is_same<typename Kokkos::View<T, P...>::memory_space,
            typename Kokkos::View<T, P...>::HostMirror::memory_space>::value &&
        std::is_same<typename Kokkos::View<T, P...>::data_type,
            typename Kokkos::View<T, P...>::HostMirror::data_type>::value)>::
        type* = 0) {
  return src;
}

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_uninit_mirror_view(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<!(
        std::is_same<typename Kokkos::View<T, P...>::memory_space,
            typename Kokkos::View<T, P...>::HostMirror::memory_space>::value &&
        std::is_same<typename Kokkos::View<T, P...>::data_type,
            typename Kokkos::View<T, P...>::HostMirror::data_type>::value)>::
        type* = 0) {
  return create_uninit_mirror(src);
}
#endif

template <typename T>
HostWrite<T>::HostWrite(LO size_in, std::string const& name)
    : write_(size_in, name)
#ifdef OMEGA_H_USE_KOKKOSCORE
      ,
      mirror_(create_uninit_mirror_view(write_.view()))
#endif
{
}

template <typename T>
HostWrite<T>::HostWrite(LO size_in, T offset, T stride, std::string const& name)
    : HostWrite<T>(Write<T>(size_in, offset, stride, name)) {}

template <typename T>
HostWrite<T>::HostWrite(Write<T> write_in)
    : write_(write_in)
#ifdef OMEGA_H_USE_KOKKOSCORE
      ,
      mirror_(create_uninit_mirror_view(write_.view()))
#endif
{
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::deep_copy(mirror_, write_.view());
#endif
}

template <typename T>
HostWrite<T>::HostWrite(std::initializer_list<T> l, std::string const& name)
    :  // an initializer_list should never have over 2 billion items...
      HostWrite<T>(static_cast<LO>(l.size()), name) {
  LO i = 0;
  for (auto v : l) operator[](i++) = v;
}

template <typename T>
Write<T> HostWrite<T>::write() const {
#ifdef OMEGA_H_USE_KOKKOSCORE
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
#ifdef OMEGA_H_USE_KOKKOSCORE
  return mirror_.data();
#else
  return write_.data();
#endif
}

template <typename T>
void HostWrite<T>::set(LO i, T value) {
#ifdef OMEGA_H_USE_KOKKOSCORE
  mirror_[i] = value;
#else
  write_[i] = value;
#endif
}

template <typename T>
T HostWrite<T>::get(LO i) const {
#ifdef OMEGA_H_USE_KOKKOSCORE
  return mirror_[i];
#else
  return write_[i];
#endif
}

template <typename T>
#ifdef __INTEL_COMPILER
HostRead<T>::HostRead() {
}
#else
HostRead<T>::HostRead() = default;
#endif

template <typename T>
HostRead<T>::HostRead(Read<T> read)
    : read_(read)
#ifdef OMEGA_H_USE_KOKKOSCORE
      ,
      mirror_(Kokkos::create_mirror_view(read.view()))
#endif
{
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::deep_copy(mirror_, read_.view());
#endif
}

template <typename T>
LO HostRead<T>::size() const {
  return read_.size();
}

#ifdef OMEGA_H_USE_KOKKOSCORE
template <typename T>
Kokkos::View<T*> Write<T>::view() const {
  return view_;
}
#endif

template <typename T>
T const* HostRead<T>::data() const {
#ifdef OMEGA_H_USE_KOKKOSCORE
  return mirror_.data();
#else
  return read_.data();
#endif
}

template <typename T>
T HostRead<T>::get(LO i) const {
  return operator[](i);
}

template <typename T>
T HostRead<T>::last() const {
  return get(size() - 1);
}

template <class T>
Write<T> deep_copy(Read<T> a) {
  Write<T> b(a.size());
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::deep_copy(b.view(), a.view());
#else
  auto f = OMEGA_H_LAMBDA(LO i) { b[i] = a[i]; };
  parallel_for(b.size(), f, "deep_copy");
#endif
  return b;
}

#define INST(T)                                                                \
  template T* nonnull(T*);                                                     \
  template T const* nonnull(T const*);                                         \
  template class NonNullPtr<T>;                                                \
  template class Write<T>;                                                     \
  template class Read<T>;                                                      \
  template class HostWrite<T>;                                                 \
  template class HostRead<T>;                                                  \
  template Write<T> deep_copy(Read<T> a);

INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
