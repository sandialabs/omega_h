#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Reals a = {1,0,2};
  bool res = is_sorted(a);
  printf("sorted: %s\n", res ? "yes" : "no");
  OMEGA_H_CHECK(res == false);
  return 0;
}
