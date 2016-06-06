#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  static Int const nx = 1;
  static Int const dim = 2;
  build_box(mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
  classify_by_angles(mesh, PI / 4);
  vtk::FullWriter writer(mesh, "out");
  writer.write();
  for (Int i = 0; i < 6; ++i) {
    mesh.add_tag(EDGE, "candidate", 1, OSH_DONT_TRANSFER,
        Read<I8>(mesh.nents(EDGE), 1));
    bool did = refine(mesh, 0.0);
    CHECK(did);
    writer.write();
  }
  }
  fini();
}
