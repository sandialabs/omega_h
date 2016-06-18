#include "omega_h.hpp"

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  OSH_CHECK(argc == 3 || argc == 4);
  osh::Mesh mesh;
  osh::binary::read(argv[1], lib.world(), &mesh);
  auto dim = mesh.dim();
  if (argc == 4) dim = atoi(argv[2]);
  osh::vtk::write_parallel(argv[argc - 1], &mesh, dim);
}
