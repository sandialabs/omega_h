#include "Omega_h_compare.hpp"
#include "access.hpp"
#include "internal.hpp"
#include "refine.hpp"

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
  mesh.add_tag<Real>(VERT, "size", 1);
  auto opts = AdaptOpts(&mesh);
  opts.min_quality_allowed = 0.47;
  do {
    Write<Real> size(mesh.nverts());
    auto coords = mesh.coords();
    auto f = LAMBDA(LO v) {
      auto x = get_vector<3>(coords, v);
      auto coarse = 0.4;
      auto fine = 0.04;
      auto radius = norm(x);
      auto diagonal = sqrt(3) - 0.5;
      auto distance = fabs(radius - 0.5) / diagonal;
      size[v] = coarse * distance + fine * (1.0 - distance);
    };
    parallel_for(mesh.nverts(), f);
    mesh.set_tag(VERT, "size", Reals(size));
    mesh.ask_lengths();
    mesh.ask_qualities();
  } while (refine_by_size(&mesh, opts));
  bool ok = check_regression("gold_corner", &mesh);
  if (!ok) return 2;
  return 0;
}
