#include <Kokkos_Core.hpp>
#include <Omega_h_library.hpp>

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  auto lib = Omega_h::Library(&argc, &argv);
  return 0;
}
