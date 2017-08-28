#ifndef OMEGA_H_ARRAY_HPP
#define OMEGA_H_ARRAY_HPP

#include <Omega_h_defines.hpp>
#include <Omega_h_kokkos.hpp>
#include <initializer_list>
#include <memory>

namespace Omega_h {

template <typename T>
class HostWrite;

template <typename T>
class Write {
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::View<T*> view_;
#else
  std::shared_ptr<T> ptr_;
  LO size_;
#endif

 public:
  OMEGA_H_INLINE Write();
#ifdef OMEGA_H_USE_KOKKOSCORE
  Write(Kokkos::View<T*> view);
#endif
  Write(LO size, std::string const& name = "");
  Write(LO size, T value, std::string const& name = "");
  Write(LO size, T offset, T stride, std::string const& name = "");
  Write(HostWrite<T> host_write);
  OMEGA_H_INLINE Write(Write<T> const& other)
      :
#ifdef OMEGA_H_USE_KOKKOSCORE
        view_(other.view_)
#else
        ptr_(other.ptr_),
        size_(other.size_)
#endif
  {
  }
  Write<T>& operator=(Write<T> const&);
  OMEGA_H_INLINE ~Write() {
#ifndef __CUDA_ARCH__
    check_release();
#endif
  }
  LO size() const;
  OMEGA_H_DEVICE T& operator[](LO i) const {
#ifdef OMEGA_H_CHECK_BOUNDS
    OMEGA_H_CHECK(0 <= i);
    OMEGA_H_CHECK(i < size());
#endif
#ifdef OMEGA_H_USE_KOKKOSCORE
    return view_(i);
#else
    return ptr_.get()[i];
#endif
  }
  OMEGA_H_INLINE T* data() const {
#ifdef OMEGA_H_USE_KOKKOSCORE
    return view_.data();
#else
    return ptr_.get();
#endif
  }
  T* nonnull_data();
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::View<T*> view() const;
#endif
  void set(LO i, T value) const;
  T get(LO i) const;
  OMEGA_H_INLINE long use_count() const {
#ifdef OMEGA_H_USE_KOKKOSCORE
    return view_.use_count();
#else
    return ptr_.use_count();
#endif
  }
  OMEGA_H_INLINE bool exists() const { return use_count() != 0; }
  std::size_t bytes() const;

 private:
  void log_allocation() const;
  void check_release() const;
};

std::size_t get_current_bytes();
std::size_t get_max_bytes();

template <typename T>
OMEGA_H_INLINE Write<T>::Write()
    :
#ifdef OMEGA_H_USE_KOKKOSCORE
      view_()
#else
      ptr_(),
      size_(0)
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
  LO size() const;
  OMEGA_H_DEVICE T const& operator[](LO i) const { return write_[i]; }
  OMEGA_H_INLINE T const* data() const { return write_.data(); }
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::View<const T*> view() const;
#endif
  T get(LO i) const;
  T first() const;
  T last() const;
  OMEGA_H_INLINE bool exists() const { return write_.exists(); }
};

class Bytes : public Read<Byte> {
 public:
  OMEGA_H_INLINE Bytes() {}
  OMEGA_H_INLINE Bytes(Read<Byte> base) : Read<Byte>(base) {}
  Bytes(Write<Byte> write);
  Bytes(LO size, Byte value, std::string const& name = "");
  Bytes(std::initializer_list<Byte> l, std::string const& name = "");
};

class LOs : public Read<LO> {
 public:
  OMEGA_H_INLINE LOs() {}
  OMEGA_H_INLINE LOs(Read<LO> base) : Read<LO>(base) {}
  LOs(Write<LO> write);
  LOs(LO size, LO value, std::string const& name = "");
  LOs(LO size, LO offset, LO stride, std::string const& name = "");
  LOs(std::initializer_list<LO> l, std::string const& name = "");
};

class GOs : public Read<GO> {
 public:
  OMEGA_H_INLINE GOs() {}
  OMEGA_H_INLINE GOs(Read<GO> base) : Read<GO>(base) {}
  GOs(Write<GO> write);
  GOs(LO size, GO value, std::string const& name = "");
  GOs(LO size, GO offset, GO stride, std::string const& name = "");
  GOs(std::initializer_list<GO> l, std::string const& name = "");
};

class Reals : public Read<Real> {
 public:
  OMEGA_H_INLINE Reals() {}
  OMEGA_H_INLINE Reals(Read<Real> base) : Read<Real>(base) {}
  Reals(Write<Real> write);
  Reals(LO size, Real value, std::string const& name = "");
  Reals(std::initializer_list<Real> l, std::string const& name = "");
};

template <typename T>
class HostRead {
  Read<T> read_;
#ifdef OMEGA_H_USE_KOKKOSCORE
  typename Kokkos::View<const T*>::HostMirror mirror_;
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
  T const* nonnull_data() const;
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
  HostWrite(LO size, std::string const& name = "");
  HostWrite(LO size, T offset, T stride, std::string const& name = "");
  HostWrite(Write<T> write);
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
  T* nonnull_data() const;
  OMEGA_H_INLINE bool exists() const { return write_.exists(); }
};

template <class T>
Write<T> deep_copy(Read<T> a);

/* begin explicit instantiation declarations */
#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template class Read<T>;                                               \
  extern template class Write<T>;                                              \
  extern template class HostRead<T>;                                           \
  extern template class HostWrite<T>;                                          \
  extern template Write<T> deep_copy(Read<T> a);
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL
/* end explicit instantiation declarations */

}  // end namespace Omega_h

#endif
