#include "omega_h.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 3);
  Omega_h::Mesh mesh;
  Omega_h::gmsh::read(argv[1], lib, &mesh);
  Omega_h::binary::write(argv[2], &mesh);
}
