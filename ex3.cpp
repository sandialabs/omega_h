#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  init(argc, argv);
  {
  auto world = Comm::world();
  Mesh mesh;
  static Int const dim = 2;
  if (world->rank() == 0) {
    build_box(mesh, 1, 1, 1, 1, 1, (dim == 3) ? 1 : 0);
    classify_by_angles(mesh, PI / 4);
  }
  vtk::FullWriter writer(mesh, "out");
  mesh.add_tag(EDGE, "candidate", 1, OSH_DONT_TRANSFER,
      Read<I8>(mesh.nedges(), 1));
  writer.write();
  refine(mesh, 0.3);
  writer.write();
  coarsen_verts(mesh, Read<I8>(mesh.nverts(), 1), 0.3, false);
  writer.write();
  }
  fini();
}
