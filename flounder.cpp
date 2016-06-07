#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 0, 1, 1, 0);
  classify_by_angles(mesh, PI / 4);
  do {
    Write<Real> size(mesh.nverts());
    auto coords = mesh.coords();
    auto f = LAMBDA(LO v) {
      auto x = get_vec<2>(coords, v);
      auto s = cos(x[0] * 8.0 * M_PI) / 4.0 + 1.0 / 2.0;
      auto d = fabs(x[1] - s);
      auto area = 1e-6 + d * 1e-4;
      auto length = sqrt(area);
      auto fudge = 1.4;
      size[v] = length * fudge;
    };
    parallel_for(mesh.nverts(), f);
    mesh.add_tag(VERT, "size", 1, OSH_DONT_TRANSFER, Reals(size));
    std::cerr << mesh.nelems() << '\n';
  } while(refine_by_size(mesh, 0.3));
  vtk::write_vtu("out.vtu", mesh, mesh.dim());
  }
  fini();
}

