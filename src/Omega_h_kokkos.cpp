#include "Omega_h_kokkos.hpp"

namespace Omega_h {

void begin_code(std::string const& name) {
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::Profiling::pushRegion(name);
#else
  (void)name;
#endif
}

void end_code() {
#ifdef OMEGA_H_USE_KOKKOSCORE
  Kokkos::Profiling::popRegion();
#endif
}

}  // namespace Omega_h
