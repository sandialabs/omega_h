#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 1, 10, 10, 10);
  Reals coords = mesh.coords();
  LOs sorted2input = hilbert::sort_coords(coords, mesh.dim());
  LOs nums = invert_permutation(sorted2input);
  mesh.add_tag<LO>(VERT, "hilbert", 1);
  mesh.set_tag<LO>(VERT, "hilbert", nums);
  std::ofstream file("out.vtu");
  vtk::write_vtu(file, mesh);
  }
  fini();
}
