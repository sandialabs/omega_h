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
  //auto coords = mesh.coords();
    auto f = LAMBDA(LO v) {
    //auto x = get_vec<2>(coords, v);
    //double s = cos(x[0] * 8.0 * PI) / 4.0 + 1.0 / 2.0;
    //double d = fabs(x[1] - s);
    //double fudge = 10;
    //size[v] = sqrt(2 * (1e-7 + d * 1e-5)) * fudge;
      size[v] = 0.3;
    };
    parallel_for(mesh.nverts(), f);
    mesh.add_tag(VERT, "size", 1, OSH_DONT_TRANSFER, Reals(size));
    writer.write();
  } while(refine_by_size(mesh, 0.3));
  }
  fini();
}

