#include "Omega_h_kokkos.hpp"

namespace Omega_h {

void begin_code(const char* name) {
#if 0
  Kokkos::Profiling::pushRegion(name);
#else
  (void)name;
#endif
}

void end_code() {
#if 0
  Kokkos::Profiling::popRegion();
#endif
}

}  // namespace Omega_h
