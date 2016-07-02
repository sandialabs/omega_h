#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  constexpr Int dim = 2;
  Mesh mesh;
  if (world->rank() == 0) {
    auto nx = 4;
    build_box(&mesh, lib, 1, 1, 0, nx, nx, (dim == 3) ? nx : 0);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  auto metric = compose_metric(identity_matrix<2,2>(),
      vector_2(1.0 / 32.0, 1.0 / 4.0));
  auto metrics = repeat_symm(mesh.nverts(), metric);
  mesh.add_tag(VERT, "metric", symm_dofs(mesh.dim()),
      OSH_METRIC, metrics);
  mesh.ask_edge_lengths();
  mesh.ask_qualities();
  vtk::FullWriter writer(&mesh, "out");
  writer.write();
  adapt(&mesh, 0.10, 0.20, 2.0 / 3.0, 4.0 / 3.0, 4, 2);
  writer.write();
}
