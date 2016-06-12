#include "internal.hpp"

using namespace osh;

static void serial_test(Mesh* mesh) {
  static Int const nx = 5;
  static Int const dim = 3;
  build_box(mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
  classify_by_angles(mesh, PI / 4);
  mesh.reorder();
  mesh.forget_globals();
  Vector<dim> lengths;
  for (Int i = 0; i < dim; ++i)
    lengths[i] = 1.0;
  lengths[dim - 1] = 0.5;
  auto metric = compose_metric(identity_matrix<dim,dim>(), lengths);
  mesh.add_tag(VERT, "metric", symm_dofs(dim), OSH_METRIC,
      repeat_symm(mesh.nverts(), metric));
  mesh.ask_edge_lengths();
  mesh.ask_qualities();
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto self = lib.self();
  Mesh mesh;
  if (world->rank() == 0) {
    serial_test(mesh);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_partition(VERTEX_BASED);
  if (mesh.dim() == 3) {
  vtk::write_parallel("tets", mesh, 3);
  }
  vtk::write_parallel("tris", mesh, 2);
  vtk::write_parallel("edges", mesh, 1);
}
