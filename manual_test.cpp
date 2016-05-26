#include <fstream>

#include "internal.hpp"

int main(int argc, char** argv) {
  init(argc, argv);
  {
  Mesh mesh;
  static Int const nx = 2;
  static Int const dim = 3;
  build_box(mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
  classify_by_angles(mesh, PI / 4);
  auto coords = mesh.coords();
  LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh.dim());
  for (Int d = 0; d <= mesh.dim(); ++d) {
    mesh.ask_owners(d);
  }
  reorder_mesh(mesh, new_verts2old_verts);
  for (Int d = 0; d <= mesh.dim(); ++d) {
    CHECK(Read<I32>(mesh.nents(d), 0) == mesh.ask_own_ranks(d));
    CHECK(LOs(mesh.nents(d), 0, 1) == mesh.ask_own_idxs(d));
  }
  mesh.forget_globals();
  auto metric = compose_metric(identity_matrix<dim,dim>(), vector_3(1,1,.5));
  mesh.add_tag(VERT, "metric", symm_dofs(dim), repeat_symm(mesh.nverts(), metric));
  mesh.ask_edge_lengths();
  mesh.ask_qualities();
  if (mesh.dim() == 3) {
  std::ofstream file("tets.vtu");
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
