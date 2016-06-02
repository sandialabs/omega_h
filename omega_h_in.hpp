#ifndef OMEGA_H_HPP
#define OMEGA_H_HPP

#cmakedefine OSH_USE_MPI
#cmakedefine OSH_USE_KOKKOS
#cmakedefine OSH_USE_OPENMP
#cmakedefine OSH_USE_CUDA
#cmakedefine OSH_USE_ZLIB
#cmakedefine OSH_CHECK_BOUNDS

#include <memory>

#ifdef OSH_USE_MPI

/* on BlueGene/Q the default install
 * defines MPICH2_CONST to empty if we
 * don't define it first, causing tons
 * of compile errors.
 *
 * in addition, the mpicxx.h header is full
 * of "const MPICH2_CONST", probably as a workaround
 * for the first mistake above.
 * as a result, properly defining MPICH2_CONST to const
 * also causes compile errors.
 * luckily, we can avoid including mpicxx.h with
 * MPICH_SKIP_MPICXX.
 */
#ifdef __bgq__
#define MPICH2_CONST const
#define MPICH_SKIP_MPICXX
#endif

/* have Clang diagnostics ignore everything
 * inside mpi.h
 */
#ifdef __clang__
#pragma clang system_header
#endif
#include <mpi.h>

#endif //OSH_USE_MPI

#ifdef OSH_USE_KOKKOS
#ifdef __clang__
#pragma clang system_header
#endif
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
#include <Kokkos_Core.hpp>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
#endif //OSH_USE_KOKKOS

#ifdef OSH_USE_KOKKOS
#define OSH_INLINE KOKKOS_INLINE_FUNCTION
#define OSH_LAMBDA KOKKOS_LAMBDA
#else
#define OSH_INLINE inline
#define OSH_LAMBDA [=]
#endif //OSH_USE_KOKKOS

#ifdef OSH_USE_CUDA
#define OSH_DEVICE __device__ inline
#else
#define OSH_DEVICE inline
#endif //OSH_USE_CUDA

//namespace osh will start here

typedef std::int8_t  I8;
typedef std::int16_t I16;
typedef std::int32_t I32;
typedef std::int64_t I64;
typedef I32          Int;
typedef I32          LO;
typedef I64          GO;
typedef double       Real;

template <typename T>
class HostWrite;

template <typename T>
class Write {
#ifdef OSH_USE_KOKKOS
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
  OSH_INLINE T& operator[](LO i) const {
#ifdef OSH_CHECK_BOUNDS
    if (i < 0)
      std::cerr << "i = " << i << '\n';
    CHECK(0 <= i);
    if (i >= size())
      std::cerr << "i = " << i << '\n';
    CHECK(i < size());
#endif
#ifdef OSH_USE_KOKKOS
    return view_(i);
#else
    return ptr_.get()[i];
#endif
  }
  T* data() const;
#ifdef OSH_USE_KOKKOS
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
  OSH_INLINE T const& operator[](LO i) const {
    return write_[i];
  }
  T const* data() const;
#ifdef OSH_USE_KOKKOS
  Kokkos::View<const T*> view() const;
#endif
  T get(LO i) const;
  T last() const;
  bool exists() const;
};

class LOs : public Read<LO> {
public:
  LOs();
  LOs(Read<LO> base);
  LOs(Write<LO> write);
  LOs(LO size, LO value);
  LOs(LO size, LO offset, LO stride);
  LOs(std::initializer_list<LO> l);
};

//namespace osh will end here

#endif
