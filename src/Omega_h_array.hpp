#ifndef OMEGA_H_ARRAY_HPP
#define OMEGA_H_ARRAY_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_fail.hpp>
#include <initializer_list>
#include <memory>
#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_kokkos.hpp>
#else
#include <string>
#endif

namespace Omega_h {

template <typename T>
T* nonnull(T* p);

template <typename T>
class HostWrite;

#ifndef OMEGA_H_USE_KOKKOSCORE
template <typename T>
struct SharedAlloc2 {
  std::string name;
  std::unique_ptr<T[]> ptr;
};
#endif

template <typename T>
class Write {
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::View<T*> view_;
#else
  std::shared_ptr<SharedAlloc2<T> > tracker_;
  LO size_;
  T* ptr_;
#endif

 public:
  OMEGA_H_INLINE Write();
#ifdef OMEGA_H_USE_KOKKOSCORE
  Write(Kokkos::View<T*> view_in);
#endif
  Write(LO size_in, std::string const& name = "");
  Write(LO size_in, T value, std::string const& name = "");
  Write(LO size_in, T offset, T stride, std::string const& name = "");
  Write(HostWrite<T> host_write);
  OMEGA_H_INLINE LO size() const {
    OMEGA_H_CHECK(exists());
#ifdef OMEGA_H_USE_KOKKOSCORE
    return static_cast<LO>(view_.size());
#else
    return size_;
#endif
  }
  OMEGA_H_DEVICE T& operator[](LO i) const {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(0 <= i);
    OMEGA_H_CHECK(i < size());
#endif
#ifdef OMEGA_H_USE_KOKKOSCORE
    return view_(i);
#else
    return ptr_[i];
#endif
  }
  OMEGA_H_INLINE T* data() const {
#ifdef OMEGA_H_USE_KOKKOSCORE
    return view_.data();
#else
    return ptr_;
#endif
  }
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::View<T*> view() const;
#endif
  void set(LO i, T value) const;
  T get(LO i) const;
  OMEGA_H_INLINE long use_count() const {
#ifdef OMEGA_H_USE_KOKKOSCORE
    return view_.use_count();
#else
    return tracker_.use_count();
#endif
  }
  OMEGA_H_INLINE bool exists() const {
#ifdef __CUDA_ARCH__
    return true;
#else
    return use_count() != 0;
#endif
  }
#ifdef OMEGA_H_USE_KOKKOSCORE
  std::string name() const;
#else
  std::string const& name() const { return tracker_->name; }
#endif
#ifndef OMEGA_H_USE_KOKKOSCORE
  void rename(std::string const& new_name);
#endif
};

template <typename T>
OMEGA_H_INLINE Write<T>::Write()
    :
#ifdef OMEGA_H_USE_KOKKOSCORE
      view_()
#else
      tracker_(),
      size_(0),
      ptr_(nullptr)
#endif
{
}

template <typename T>
class Read {
  Write<T> write_;

 public:
  OMEGA_H_INLINE Read() {}
  Read(Write<T> write);
  Read(LO size, T value, std::string const& name = "");
  Read(LO size, T offset, T stride, std::string const& name = "");
  Read(std::initializer_list<T> l, std::string const& name = "");
  OMEGA_H_INLINE LO size() const { return write_.size(); }
  OMEGA_H_DEVICE T const& operator[](LO i) const { return write_[i]; }
  OMEGA_H_INLINE T const* data() const { return write_.data(); }
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::View<const T*> view() const;
#endif
  T get(LO i) const;
  T first() const;
  T last() const;
  OMEGA_H_INLINE bool exists() const { return write_.exists(); }
  std::string name() const { return write_.name(); }
};

template <typename T>
Read<T> read(Write<T> a) {
  return Read<T>(a);
}

class Bytes : public Read<Byte> {
 public:
  OMEGA_H_INLINE Bytes() {}
  OMEGA_H_INLINE Bytes(Read<Byte> base) : Read<Byte>(base) {}
  Bytes(Write<Byte> write);
  Bytes(LO size_in, Byte value, std::string const& name = "");
  Bytes(std::initializer_list<Byte> l, std::string const& name = "");
};

class LOs : public Read<LO> {
 public:
  OMEGA_H_INLINE LOs() {}
  OMEGA_H_INLINE LOs(Read<LO> base) : Read<LO>(base) {}
  LOs(Write<LO> write);
  LOs(LO size_in, LO value, std::string const& name = "");
  LOs(LO size_in, LO offset, LO stride, std::string const& name = "");
  LOs(std::initializer_list<LO> l, std::string const& name = "");
};

class GOs : public Read<GO> {
 public:
  OMEGA_H_INLINE GOs() {}
  OMEGA_H_INLINE GOs(Read<GO> base) : Read<GO>(base) {}
  GOs(Write<GO> write);
  GOs(LO size_in, GO value, std::string const& name = "");
  GOs(LO size_in, GO offset, GO stride, std::string const& name = "");
  GOs(std::initializer_list<GO> l, std::string const& name = "");
};

class Reals : public Read<Real> {
 public:
  OMEGA_H_INLINE Reals() {}
  OMEGA_H_INLINE Reals(Read<Real> base) : Read<Real>(base) {}
  Reals(Write<Real> write);
  Reals(LO size_in, Real value, std::string const& name = "");
  Reals(std::initializer_list<Real> l, std::string const& name = "");
};

template <typename T>
class HostRead {
  Read<T> read_;
#ifdef OMEGA_H_USE_KOKKOSCORE
  typename Kokkos::View<const T*, Kokkos::HostSpace> mirror_;
#endif
 public:
  HostRead();
  HostRead(Read<T> read);
  LO size() const;
  inline T const& operator[](LO i) const {
#ifdef OMEGA_H_USE_KOKKOSCORE
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(0 <= i);
    OMEGA_H_CHECK(i < size());
#endif
    return mirror_(i);
#else
    return read_[i];
#endif
  }
  T const* data() const;
  T get(LO i) const;
  T last() const;
};

template <typename T>
class HostWrite {
  Write<T> write_;
#ifdef OMEGA_H_USE_KOKKOSCORE
  typename Kokkos::View<T*>::HostMirror mirror_;
#endif
 public:
  HostWrite();
  HostWrite(LO size_in, std::string const& name = "");
  HostWrite(LO size_in, T offset, T stride, std::string const& name = "");
  HostWrite(Write<T> write_in);
  HostWrite(std::initializer_list<T> l, std::string const& name = "");
  Write<T> write() const;
  LO size() const;
  inline T& operator[](LO i) const {
#ifdef OMEGA_H_USE_KOKKOSCORE
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(0 <= i);
    OMEGA_H_CHECK(i < size());
#endif
    return mirror_(i);
#else
    return write_[i];
#endif
  }
  T* data() const;
  OMEGA_H_INLINE bool exists() const { return write_.exists(); }
  void set(LO i, T value);
  T get(LO i) const;
};

template <class T>
void copy_into(Read<T> a, Write<T> b);
template <class T>
Write<T> deep_copy(Read<T> a, std::string const& name = "");

/* begin explicit instantiation declarations */
#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template T* nonnull(T*);                                              \
  extern template T const* nonnull(T const*);                                  \
  extern template class Read<T>;                                               \
  extern template class Write<T>;                                              \
  extern template class HostRead<T>;                                           \
  extern template class HostWrite<T>;                                          \
  extern template void copy_into(Read<T> a, Write<T> b); \
  extern template Write<T> deep_copy(Read<T> a, std::string const&);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL
/* end explicit instantiation declarations */

}  // end namespace Omega_h

#endif
