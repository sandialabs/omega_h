#include "Omega_h_for.hpp"
#include "Omega_h_library.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  const int size = 5;
  Write<Real> out(size);
  auto f = OMEGA_H_LAMBDA(LO i) {
    out[i] = i;
  };
  parallel_for(size, f, "foo");
  HostRead<Real> hr(out);
  for(int i=0; i<size; i++)
    assert(hr[i]==i);
  return 0;
}
