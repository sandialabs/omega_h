#include <Omega_h_build.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_timer.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>

#include <iostream>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto mesh = build_box(world, 1., 1., 0.5, 8, 8, 4);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto metrics = get_implied_metrics(&mesh);
  mesh.add_tag(VERT, "metric", symm_ncomps(mesh.dim()), metrics);
  auto target_metric = compose_metric(
      identity_matrix<3, 3>(), vector_3(1.0 / 64.0, 1.0 / 4.0, 1.0 / 8.0));
  auto target_metrics = repeat_symm(mesh.nverts(), target_metric);
  mesh.add_tag(VERT, "target_metric", symm_ncomps(mesh.dim()), target_metrics);
  mesh.ask_lengths();
  mesh.ask_qualities();
  auto opts = AdaptOpts(&mesh);
  opts.verbosity = EXTRA_STATS;
  Now t0 = now();
  while (approach_metric(&mesh, opts)) adapt(&mesh, opts);
  Now t1 = now();
  std::cout << "anisotropic approach took " << (t1 - t0) << " seconds\n";
  bool ok = check_regression("gold_aniso", &mesh);
  if (!ok) return 2;
  return 0;
}
