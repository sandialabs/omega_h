#include "coarsen.hpp"
#include "Omega_h.hpp"

#include <cmath>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh;
  if (world->rank() == 0) {
    build_box(&mesh, lib, 1, 1, 1, 4, 4, 4);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.add_tag<Real>(VERT, "size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT);
  mesh.set_tag(VERT, "size", Reals(mesh.nverts(), 1.0));
  auto opts = AdaptOpts(&mesh);
  opts.min_quality_allowed = 0.47;
  while (coarsen_by_size(&mesh, opts))
    ;
  bool ok = check_regression("gold_coarsen", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}
