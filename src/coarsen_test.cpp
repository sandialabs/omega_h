#include "Omega_h.hpp"
#include "Omega_h_coarsen.hpp"
#include "Omega_h_compare.hpp"

#include <cmath>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh(&lib);
  if (world->rank() == 0) {
    build_box(&mesh, 1, 1, 1, 4, 4, 4);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(VERT, "metric", Reals(mesh.nverts(), 1.0));
  auto opts = AdaptOpts(&mesh);
  opts.min_quality_allowed = 0.47;
  while (coarsen_by_size(&mesh, opts)) {
  }
  bool ok = check_regression("gold_coarsen", &mesh);
  if (!ok) return 2;
  return 0;
}
