#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>

#include <cmath>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, 1., 1., 1., 4, 4, 4);
  auto opts = AdaptOpts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(0.3)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(0.6)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.0)));
  while (coarsen_by_size(&mesh, opts))
    ;
  mesh.ask_qualities();
  bool ok = check_regression("gold_coarsen", &mesh);
  if (!ok) return 2;
  return 0;
}
