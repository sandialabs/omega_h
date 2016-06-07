#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 0, 1, 1, 0);
  classify_by_angles(mesh, PI / 4);
  vtk::FullWriter writer(mesh, "out");
  do {
    Write<Real> size(mesh.nverts());
    auto coords = mesh.coords();
    auto f = LAMBDA(LO v) {
      auto x = get_vec<2>(coords, v);
      auto area = 1e-5 + (1e-3) * x[0];
      auto length = sqrt(area);
      size[v] = length;
    };
    parallel_for(mesh.nverts(), f);
    mesh.add_tag(VERT, "size", 1, OSH_DONT_TRANSFER, Reals(size));
    std::cerr << mesh.nelems() << '\n';
    writer.write();
  } while(refine_by_size(mesh, 0.3));
  }
  fini();
}

