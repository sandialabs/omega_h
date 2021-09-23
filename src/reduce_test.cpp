#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  const int size = 5;
  Write<Real> out(size);

  auto transform = OMEGA_H_LAMBDA(LO i)->bool {
    return true;
  };

  bool res;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<>(0, 1), 
    KOKKOS_LAMBDA(const int& i, bool& update) {
      update = transform(i);
  }, Kokkos::LAnd<bool>(res) );
  assert(res==true);
  return 0;
}
