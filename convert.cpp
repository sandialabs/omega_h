#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  CHECK(argc == 3);
  Mesh mesh;
  {
  std::ifstream in(argv[1]);
  gmsh::read(in, mesh);
  }
  vtk::write_vtu(argv[2], mesh, mesh.dim());
  }
  fini();
}
