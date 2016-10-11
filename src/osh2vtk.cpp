#include "Omega_h.hpp"

#include <cstdlib>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 3 || argc == 4);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], lib.world(), &mesh);
  auto dim = mesh.dim();
  if (argc == 4) dim = atoi(argv[2]);
  Omega_h::vtk::write_parallel(argv[argc - 1], &mesh, dim);
}
