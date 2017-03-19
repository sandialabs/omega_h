#include "Omega_h.hpp"
#include "Omega_h_math.hpp"
#include "Omega_h_timer.hpp"

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
void run_case(Mesh* mesh, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = find_implied_metric(mesh);
  mesh->add_tag(VERT, "metric", symm_dofs(dim),
      implied_metrics);
  mesh->add_tag<Real>(
      VERT, "target_metric", symm_dofs(dim));
  set_target_metric<dim>(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->ask_lengths();
  mesh->ask_qualities();
  vtk::FullWriter writer;
  if (vtk_path) {
    writer = vtk::FullWriter(mesh, vtk_path);
    writer.write();
  }
  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  Now t0 = now();
  while (approach_size_field(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
    if (vtk_path) writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  char const* in_path = nullptr;
  char const* out_path = nullptr;
  char const* vtk_path = nullptr;
  bool should_help = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string("-a") == argv[i]) {
      if (i == argc - 1) {
        std::cout << "-a takes an argument\n";
        should_help = true;
      } else {
        vtk_path = argv[++i];
      }
    } else if (!in_path) {
      in_path = argv[i];
    } else if (!out_path) {
      out_path = argv[i];
    } else {
      std::cout << "unknown argument " << argv[i] << '\n';
      should_help = true;
    }
  }
  if (!in_path || !out_path) should_help = true;
  if (should_help) {
    std::cout << "usage: " << argv[0]
              << " [options] input.mesh[b] output.mesh[b]\n";
    std::cout
        << "options: -a vtk_path                debug output for adaptivity\n";
    return -1;
  }
  Mesh mesh(&lib);
  meshb::read(&mesh, in_path);
  if (mesh.dim() == 2) run_case<2>(&mesh, vtk_path);
  if (mesh.dim() == 3) run_case<3>(&mesh, vtk_path);
  meshb::write(&mesh, out_path, 2);
  return 0;
}
