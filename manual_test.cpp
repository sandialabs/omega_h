#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 1, 2, 2, 2);
  std::ofstream file("out.vtu");
  vtk::write_vtu(file, mesh);
  }
  fini();
}
