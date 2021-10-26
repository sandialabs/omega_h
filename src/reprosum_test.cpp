#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Reals a = {0,1,2};

  Real res = repro_sum(a);
  printf("result %f\n", res);
  OMEGA_H_CHECK(are_close(res,3.0));
  return 0;
}
