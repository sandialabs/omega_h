#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  CHECK(argc == 3);
  Mesh mesh;
  gmsh::read(argv[1], mesh);
  vtk::write_vtu(argv[2], mesh, mesh->dim());
}
