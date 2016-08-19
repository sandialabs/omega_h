#include "access.hpp"
#include "adapt.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "refine.hpp"

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
  mesh.add_tag<Real>(VERT, "size", 1, OMEGA_H_LINEAR_INTERP, OMEGA_H_DO_OUTPUT);
  do {
    Write<Real> size(mesh.nverts());
    auto coords = mesh.coords();
    auto f = LAMBDA(LO v) {
      auto x = get_vector<3>(coords, v);
      auto coarse = 0.4;
      auto fine = 0.06;
      auto radius = norm(x);
      auto diagonal = sqrt(3) - 0.5;
      auto distance = fabs(radius - 0.5) / diagonal;
      size[v] = coarse * distance + fine * (1.0 - distance);
    };
    parallel_for(mesh.nverts(), f);
    mesh.set_tag(VERT, "size", Reals(size));
    mesh.ask_lengths();
    mesh.ask_qualities();
  } while (refine_by_size(&mesh, 4.0 / 3.0, 0.47, false));
  bool ok = check_regression("gold_corner", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}
