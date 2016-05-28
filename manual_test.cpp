#include <fstream>

#include "internal.hpp"

static void serial_test(Mesh& mesh) {
  static Int const nx = 1;
  static Int const dim = 2;
  build_box(mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
  classify_by_angles(mesh, PI / 4);
  auto coords = mesh.coords();
  LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh.dim());
  reorder_mesh(mesh, new_verts2old_verts);
  mesh.forget_globals();
  Vector<dim> lengths;
  for (Int i = 0; i < dim; ++i)
    lengths[i] = 1.0;
  lengths[dim - 1] = 0.5;
  auto metric = compose_metric(identity_matrix<dim,dim>(), lengths);
  mesh.add_tag(VERT, "metric", symm_dofs(dim), repeat_symm(mesh.nverts(), metric));
  mesh.ask_edge_lengths();
  mesh.ask_qualities();
  if (mesh.dim() == 3) {
  std::ofstream file("tets.vtu");
  vtk::write_vtu(file, mesh, 3);
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

int main(int argc, char** argv) {
  init(argc, argv);
  {
  auto world = Comm::world();
  auto self = Comm::self();
  Mesh mesh;
  if (world->rank() == 0) {
    serial_test(mesh);
  } else {
    mesh.set_comm(self);
  }
  bcast_mesh(mesh, world, world->rank() == 0);
  mesh.set_comm(world);
  if (world->rank() == 0) {
    migrate_mesh(mesh, Remotes(Read<I32>({0}), LOs({0})));
  } else {
    migrate_mesh(mesh, Remotes(Read<I32>({0}), LOs({1})));
  }
  }
  fini();
}
