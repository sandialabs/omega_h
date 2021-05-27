#include <Kokkos_Core.hpp>
#include <Omega_h_library.hpp>

using namespace Omega_h;

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  auto lib = Library(&argc, &argv);
  return 0;
}
