#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 3);
  LOs rv2v = mesh.ask_adj(3,0).ab2b;
  CHECK(rv2v.size() == 4 * 6);
  for (Int i = 0; i < 6; ++i) {
    for (Int j = 0; j < 4; ++j)
      std::cerr << ' ' << rv2v[i * 4 + j];
    std::cerr << '\n';
  }
  std::ofstream file("out.vtu");
  vtk::write_vtu(file, mesh);
  }
  fini();
}
