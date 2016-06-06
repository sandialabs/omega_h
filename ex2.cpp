#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  init(argc, argv);
  {
  auto world = Comm::world();
  Mesh mesh;
  static Int const dim = 2;
  if (world->rank() == 0) {
    build_box(mesh, 2, 1, 1, 2, 1, (dim == 3) ? 1 : 0);
    classify_by_angles(mesh, PI / 4);
  }
//mesh.set_comm(world);
//mesh.balance();
  vtk::FullWriter writer(mesh, "out");
  writer.write();
  for (Int i = 0; i < 1; ++i) {
    mesh.add_tag(EDGE, "candidate", 1, OSH_DONT_TRANSFER,
        Read<I8>(mesh.nents(EDGE), 1));
    bool did = refine(mesh, 0.0);
    CHECK(did);
    writer.write();
  }
  }
  fini();
}
