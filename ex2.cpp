#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  static Int const nx = 1;
  static Int const dim = 2;
  build_box(mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
  classify_by_angles(mesh, PI / 4);
  mesh.add_tag(EDGE, "candidate", 1, OSH_DONT_TRANSFER,
      Read<I8>(mesh.nents(EDGE), 1));
  bool did = refine(mesh, 0.0);
  CHECK(did);
  if (mesh.dim() == 3) {
  vtk::write_parallel("tets", mesh, 3);
  }
  vtk::write_parallel("tris", mesh, 2);
  vtk::write_parallel("edges", mesh, 1);
  }
  fini();
}
