#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 1, 1, 1, 1);
  Reals coords = mesh.coords();
  LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh.dim());
  reorder_mesh(mesh, new_verts2old_verts);
//{
//std::ofstream file("tets.vtu");
//vtk::write_vtu(file, mesh, 3);
//}
//{
//std::ofstream file("tris.vtu");
//vtk::write_vtu(file, mesh, 2);
//}
//{
//std::ofstream file("edges.vtu");
//vtk::write_vtu(file, mesh, 1);
//}
  }
  fini();
}
