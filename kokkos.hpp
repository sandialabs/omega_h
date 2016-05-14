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
