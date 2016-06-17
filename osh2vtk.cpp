#include "omega_h.hpp"

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  OSH_CHECK(argc == 3);
  osh::Mesh mesh;
  osh::binary::read(argv[1], lib.world(), &mesh);
  osh::vtk::write_parallel(argv[2], &mesh, mesh.dim());
}
