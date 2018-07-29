#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_library.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_refine.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, 1., 1., 1., 4, 4, 4);
  mesh.add_tag<Real>(VERT, "metric", 1);
  auto opts = AdaptOpts(&mesh);
  opts.min_quality_allowed = 0.47;
  do {
    Write<Real> metrics_w(mesh.nverts());
    auto coords = mesh.coords();
    auto f = OMEGA_H_LAMBDA(LO v) {
      auto x = get_vector<3>(coords, v);
      auto coarse = 0.4;
      auto fine = 0.04;
      auto radius = norm(x);
      auto diagonal = std::sqrt(3) - 0.5;
      auto distance = std::abs(radius - 0.5) / diagonal;
      auto h = coarse * distance + fine * (1.0 - distance);
      metrics_w[v] = metric_eigenvalue_from_length(h);
    };
    parallel_for(mesh.nverts(), f);
    mesh.set_tag(VERT, "metric", Reals(metrics_w));
    mesh.ask_lengths();
    mesh.ask_qualities();
  } while (refine_by_size(&mesh, opts));
  bool ok = check_regression("gold_corner", &mesh);
  if (!ok) return 2;
  return 0;
}
