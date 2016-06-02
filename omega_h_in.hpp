#ifndef OMEGA_H_HPP
#define OMEGA_H_HPP

#cmakedefine OSH_USE_MPI
#cmakedefine OSH_USE_KOKKOS
#cmakedefine OSH_USE_OPENMP
#cmakedefine OSH_USE_CUDA
#cmakedefine OSH_USE_ZLIB
#cmakedefine OSH_CHECK_BOUNDS

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

#endif
