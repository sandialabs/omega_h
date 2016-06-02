#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  static Int const nx = 1;
  static Int const dim = 2;
  build_box(mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
  classify_by_angles(mesh, PI / 4);
  if (mesh.dim() == 3) {
  vtk::write_parallel("tets", mesh, 3);
  }
  vtk::write_parallel("tris", mesh, 2);
  vtk::write_parallel("edges", mesh, 1);
  }
  fini();
}
