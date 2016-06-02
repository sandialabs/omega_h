#ifndef OMEGA_H_HPP
#define OMEGA_H_HPP

#cmakedefine OSH_USE_MPI
#cmakedefine OSH_USE_KOKKOS
#cmakedefine OSH_USE_OPENMP
#cmakedefine OSH_USE_CUDA
#cmakedefine OSH_USE_ZLIB
#cmakedefine OSH_CHECK_BOUNDS

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
#endif

#ifdef OSH_USE_KOKKOS
#define OSH_INLINE KOKKOS_INLINE_FUNCTION
#define OSH_LAMBDA KOKKOS_LAMBDA
#else
#define OSH_INLINE inline
#define OSH_LAMBDA [=]
#endif

#ifdef OSH_USE_CUDA
#define OSH_DEVICE __device__ inline
#else
#define OSH_DEVICE inline
#endif

#endif
