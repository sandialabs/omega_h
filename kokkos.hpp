#ifdef USE_KOKKOS
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

#ifdef USE_KOKKOS
#define INLINE KOKKOS_INLINE_FUNCTION
#define LAMBDA KOKKOS_LAMBDA
#else
#define INLINE inline
#define LAMBDA [=]
#endif

/* in some cases involving __constant__ data,
 * we can't have an inline function marked as
 * __host__ at all (which KOKKOS_INLINE_FUNCTION does),
 * so we need a separate macro for things that
 * execute only on the device */

#ifdef USE_CUDA
#define DEVICE __device__ inline
#else
#define DEVICE inline
#endif
