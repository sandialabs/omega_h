#include "omega_h.hpp"

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 3);
  osh::Mesh mesh;
  osh::gmsh::read(argv[1], lib, &mesh);
  osh::binary::write(argv[2], &mesh);
}
