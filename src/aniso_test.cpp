#include "omega_h.hpp"
#include "omega_h_math.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh;
  if (world->rank() == 0) {
    build_box(&mesh, lib, 1, 1, .5, 8, 8, 4);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OSH_GHOSTED);
  auto metrics = find_identity_metric(&mesh);
  mesh.add_tag(VERT, "metric", symm_dofs(mesh.dim()), OSH_METRIC, metrics);
  auto target_metric = compose_metric(
      identity_matrix<3, 3>(), vector_3(1.0 / 64.0, 1.0 / 4.0, 1.0 / 8.0));
  auto target_metrics = repeat_symm(mesh.nverts(), target_metric);
  mesh.add_tag(VERT, "target_metric", symm_dofs(mesh.dim()), OSH_METRIC,
               target_metrics);
  mesh.ask_lengths();
  mesh.ask_qualities();
  vtk::FullWriter writer(&mesh, "out");
  writer.write();
  while (approach_metric(&mesh, 0.20)) {
    adapt(&mesh, 0.20, 0.30, 2.0 / 3.0, 4.0 / 3.0, 4, 2);
    writer.write();
  }
}
