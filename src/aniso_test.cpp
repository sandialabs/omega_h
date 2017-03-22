#include "Omega_h.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_timer.hpp"

#include <iostream>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh(&lib);
  if (world->rank() == 0) {
    build_box(&mesh, 1, 1, .5, 8, 8, 4);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto metrics = find_implied_metric(&mesh);
  mesh.add_tag(VERT, "metric", symm_dofs(mesh.dim()), metrics);
  auto target_metric = compose_metric(
      identity_matrix<3, 3>(), vector_3(1.0 / 64.0, 1.0 / 4.0, 1.0 / 8.0));
  auto target_metrics = repeat_symm(mesh.nverts(), target_metric);
  mesh.add_tag(VERT, "target_metric", symm_dofs(mesh.dim()), target_metrics);
  mesh.ask_lengths();
  mesh.ask_qualities();
  auto opts = AdaptOpts(&mesh);
  Now t0 = now();
  while (approach_metric(&mesh, opts)) adapt(&mesh, opts);
  Now t1 = now();
  std::cout << "anisotropic approach took " << (t1 - t0) << " seconds\n";
  bool ok = check_regression("gold_aniso", &mesh);
  if (!ok) return 2;
  return 0;
}
