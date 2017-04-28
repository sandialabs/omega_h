#include "Omega_h_array.hpp"

#include <cstring>
#include <sstream>

#include "Omega_h_control.hpp"
#include "Omega_h_functors.hpp"
#include "Omega_h_internal.hpp"

namespace Omega_h {

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

#ifdef OMEGA_H_USE_KOKKOS
template <typename T>
Write<T>::Write(Kokkos::View<T*> view) : view_(view) {
  log_allocation();
}
#endif

template <typename T>
Write<T>::Write(LO size)
    :
#ifdef OMEGA_H_USE_KOKKOS
      view_(Kokkos::ViewAllocateWithoutInitializing("omega_h"),
          static_cast<std::size_t>(size))
#else
      ptr_(new T[size], std::default_delete<T[]>()),
      size_(size)
#endif
{
  log_allocation();
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
#ifdef OMEGA_H_USE_KOKKOS
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
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(LO size, T value) : Write<T>(size) {
  fill(*this, value);
}

template <typename T>
void fill_linear(Write<T> a, T offset, T stride) {
  auto f = OMEGA_H_LAMBDA(LO i) { a[i] = offset + (stride * static_cast<T>(i)); };
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
#ifdef OMEGA_H_USE_KOKKOS
  return static_cast<LO>(view_.size());
#else
  return size_;
#endif
}

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

Reals::Reals() : Read<Real>() {}

Reals::Reals(Write<Real> write) : Read<Real>(write) {}

Reals::Reals(LO size, Real value) : Read<Real>(size, value) {}

Reals::Reals(std::initializer_list<Real> l) : Read<Real>(l) {}

LOs::LOs(Write<LO> write) : Read<LO>(write) {}

LOs::LOs(LO size, LO value) : Read<LO>(size, value) {}

LOs::LOs(LO size, LO offset, LO stride) : Read<LO>(size, offset, stride) {}

LOs::LOs(std::initializer_list<LO> l) : Read<LO>(l) {}

template <typename T>
Read<T>::Read(Write<T> write)
    : write_(write)
#ifdef OMEGA_H_USE_KOKKOS
      ,
      access_view_(write.view())
#endif
{
}

template <typename T>
Read<T>::Read(LO size, T value) : Read<T>(Write<T>(size, value)) {}

template <typename T>
Read<T>::Read(LO size, T offset, T stride)
    : Read<T>(Write<T>(size, offset, stride)) {}

template <typename T>
Read<T>::Read(std::initializer_list<T> l) : Read<T>(HostWrite<T>(l).write()) {}

template <typename T>
LO Read<T>::size() const {
  return write_.size();
}

#ifdef OMEGA_H_USE_KOKKOS
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

template <typename T>
HostWrite<T>::HostWrite(LO size)
    : write_(size)
#ifdef OMEGA_H_USE_KOKKOS
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
#ifdef OMEGA_H_USE_KOKKOS
      ,
      mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
#ifdef OMEGA_H_USE_KOKKOS
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
#ifdef OMEGA_H_USE_KOKKOS
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
#ifdef OMEGA_H_USE_KOKKOS
  return mirror_.data();
#else
  return write_.data();
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
T* HostWrite<T>::nonnull_data() const {
  return NonNullPtr<T>::get(data());
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
#ifdef OMEGA_H_USE_KOKKOS
      ,
      mirror_(Kokkos::create_mirror_view(read.view()))
#endif
{
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::deep_copy(mirror_, read_.view());
#endif
}

template <typename T>
LO HostRead<T>::size() const {
  return read_.size();
}

#ifdef OMEGA_H_USE_KOKKOS
template <typename T>
Kokkos::View<T*> Write<T>::view() const {
  return view_;
}
#endif

template <typename T>
T const* HostRead<T>::data() const {
#ifdef OMEGA_H_USE_KOKKOS
  return mirror_.data();
#else
  return read_.data();
#endif
}

template <typename T>
T const* HostRead<T>::nonnull_data() const {
  return NonNullPtr<T>::get(data());
}

template <typename T>
T HostRead<T>::last() const {
  return operator[](size() - 1);
}

template <class T>
Write<T> deep_copy(Read<T> a) {
  Write<T> b(a.size());
#ifdef OMEGA_H_USE_KOKKOS
  Kokkos::deep_copy(b.view(), a.view());
#else
  auto f = OMEGA_H_LAMBDA(LO i) { b[i] = a[i]; };
  parallel_for(b.size(), f);
#endif
  return b;
}

#define INST(T)                                                                \
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
