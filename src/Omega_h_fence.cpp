#include <Omega_h_fence.hpp>

#ifdef OMEGA_H_USE_KOKKOSCORE
#include <Omega_h_Kokkos.hpp>
#endif

namespace Omega_h {

void fence() {
#if defined(OMEGA_H_USE_KOKKOSCORE)
  Kokkos::fence();
#elif defined(OMEGA_H_USE_CUDA)
  auto const err = cudaDeviceSynchronize();
  OMEGA_H_CHECK(err == cudaSuccess);
#endif
}

}  // namespace Omega_h
