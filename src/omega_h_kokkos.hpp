#ifndef OMEGA_H_KOKKOS_HPP
#define OMEGA_H_KOKKOS_HPP

#ifdef OSH_USE_KOKKOS

OSH_SYSTEM_HEADER

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <Kokkos_Core.hpp>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#endif  // OSH_USE_KOKKOS

#ifdef OSH_USE_KOKKOS
#define OMEGA_H_INLINE KOKKOS_INLINE_FUNCTION
#else
#define OMEGA_H_INLINE inline
#endif  // OSH_USE_KOKKOS

#ifdef OSH_USE_CUDA
#define OSH_DEVICE __device__ inline
#define OSH_LAMBDA [=] __device__
#else
#define OSH_DEVICE inline
#define OSH_LAMBDA [=]
#endif  // OSH_USE_CUDA

#endif
