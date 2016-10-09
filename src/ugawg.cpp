#include "Omega_h.hpp"
#include "Omega_h_math.hpp"
#include "timer.hpp"
#include "loop.hpp"

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 2;

static void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_dofs(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.1;
    h[dim - 1] = 0.001 + 0.198 * fabs(z - 0.5);
    auto m = diagonal(metric_eigenvalues(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh;
  if (world->rank() == 0) {
    build_box(&mesh, lib, 1, 1, 1, 8, 8, 8 * (dim - 2));
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = find_implied_metric(&mesh);
  mesh.add_tag(VERT, "metric", symm_dofs(dim), OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT, implied_metrics);
  mesh.add_tag<Real>(VERT, "target_metric", symm_dofs(mesh.dim()), OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT);
  set_target_metric(&mesh);
  mesh.ask_lengths();
  mesh.ask_qualities();
  vtk::FullWriter writer(&mesh, "debug");
  writer.write();
  auto opts = AdaptOpts();
  opts.min_quality_allowed = 0.40;
  opts.min_quality_desired = 0.50;
  opts.verbosity = EXTRA_STATS;
  Now t0 = now();
  while (approach_size_field(&mesh, opts)) {
    adapt(&mesh, opts);
    if (mesh.has_tag(VERT, "target_metric")) {
      set_target_metric(&mesh);
    }
    writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
  return 0;
}

