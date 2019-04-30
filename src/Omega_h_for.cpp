#include <Omega_h_for.hpp>

namespace Omega_h {
#if defined(OMEGA_H_USE_CUDA) && (!defined(OMEGA_H_USE_KOKKOS))
int block_size_cuda = 256;
#endif
}  // namespace Omega_h
