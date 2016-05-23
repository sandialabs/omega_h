#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 1, 4, 4, 4);
  classify_by_angles(mesh, PI / 4);
//LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh.dim());
//reorder_mesh(mesh, new_verts2old_verts);
//mesh.forget_globals();
  {
  std::ofstream file("out.vtu");
  vtk::write_vtu(file, mesh);
  }
  {
  std::ofstream file("tris.vtu");
  vtk::write_vtu(file, mesh, 2);
  }
  {
  std::ofstream file("edges.vtu");
  vtk::write_vtu(file, mesh, 1);
  }
  }
  fini();
}
