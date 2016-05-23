#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  build_box(mesh, 1, 1, 1, 2, 2, 2);
  classify_by_angles(mesh, PI / 4);
  auto coords = mesh.coords();
  LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh.dim());
  reorder_mesh(mesh, new_verts2old_verts);
  mesh.forget_globals();
  auto metric = compose_metric(identity_matrix<3,3>(), vector_3(.5,.5,.5));
  mesh.add_tag(VERT, "metric", symm_dofs(3), repeat_symm(mesh.nverts(), metric));
  mesh.ask_edge_lengths();
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
