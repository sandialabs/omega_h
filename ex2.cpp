#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh;
  static Int const dim = 2;
  if (world->rank() == 0) {
    build_box(&mesh, 1, 1, 1, 1, 1, (dim == 3) ? 1 : 0);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  vtk::FullWriter writer(&mesh, "out");
  mesh.add_tag<Real>(VERT, "size", 1, OSH_LINEAR_INTERP);
  Write<Real> size(mesh.nverts());
  auto f = LAMBDA(LO v) {
    size[v] = 0.2;
  };
  parallel_for(mesh.nverts(), f);
  mesh.set_tag(VERT, "size", Reals(size));
  mesh.ask_edge_lengths();
  do {
    std::cerr << mesh.nglobal_ents(mesh.dim()) << " elements\n";
    std::cerr << mesh.nglobal_ents(EDGE) << " edges\n";
    writer.write();
  } while(refine_by_size(&mesh, 0.3));
}
