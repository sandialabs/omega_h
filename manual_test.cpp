#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 1, 1, 1, 0);
  Reals coords = mesh.coords();
  classify_by_angles(mesh, PI / 4);
  std::ofstream file("out.vtu");
  vtk::write_vtu(file, mesh);
  }
  fini();
}
