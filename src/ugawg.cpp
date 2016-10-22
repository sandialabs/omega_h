#include "Omega_h.hpp"
#include "Omega_h_math.hpp"
#include "loop.hpp"
#include "timer.hpp"

#include <iostream>

using namespace Omega_h;

template <Int dim>
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

template <Int dim>
void run_case(Mesh* mesh) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = find_implied_metric(mesh);
  mesh->add_tag(VERT, "metric", symm_dofs(dim), OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT, implied_metrics);
  mesh->add_tag<Real>(
      VERT, "target_metric", symm_dofs(dim), OMEGA_H_METRIC, OMEGA_H_DO_OUTPUT);
  set_target_metric<dim>(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->ask_lengths();
  mesh->ask_qualities();
  vtk::FullWriter writer(mesh, "debug");
  writer.write();
  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  Now t0 = now();
  while (approach_size_field(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) {
      set_target_metric<dim>(mesh);
    }
    writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  if (argc != 2) {
    std::cout << "usage: " << argv[0] << " input.osh\n";
    return -1;
  }
  Mesh mesh(&lib);
  binary::read(argv[1], lib.world(), &mesh);
  if (mesh.dim() == 2) run_case<2>(&mesh);
  if (mesh.dim() == 3) run_case<3>(&mesh);
  return 0;
}
