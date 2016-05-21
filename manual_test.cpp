#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 0, 5, 13, 0);
  std::ofstream file("out.vtu");
  vtk::write_vtu(file, mesh, 2);
  }
  fini();
}
