#include "Omega_h.hpp"
#include "Omega_h_math.hpp"
#include "Omega_h_egads.hpp"
#include "loop.hpp"
#include "timer.hpp"

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 3;

void run_case(Mesh* mesh, Egads* eg, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_size = find_implied_size(mesh);
  mesh->add_tag(VERT, "size", 1, OMEGA_H_SIZE,
      OMEGA_H_DO_OUTPUT, implied_size);
  auto target_size = multiply_each_by(0.5, implied_size);
  mesh->add_tag<Real>(
      VERT, "target_size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT,
      target_size);
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
  opts.egads_model = eg;
  Now t0 = now();
  if (vtk_path) writer.write();
  while (approach_size_field(mesh, opts)) {
    if (vtk_path) writer.write();
    adapt(mesh, opts);
    if (vtk_path) writer.write();
    break;
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  char const* in_path = nullptr;
  char const* egads_path = nullptr;
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
    } else if (ends_with(argv[i], ".egads")) {
      egads_path = argv[i];
    } else if (!in_path) {
      in_path = argv[i];
    } else if (!out_path) {
      out_path = argv[i];
    } else {
      std::cout << "unknown argument " << argv[i] << '\n';
      should_help = true;
    }
  }
  if (!in_path || !out_path || !egads_path) should_help = true;
  if (should_help) {
    std::cout << "usage: " << argv[0]
              << " [options] input.mesh[b] [input.egads] output.mesh[b]\n";
    std::cout
        << "options: -a vtk_path                debug output for adaptivity\n";
    return -1;
  }
  Mesh mesh(&lib);
  meshb::read(&mesh, in_path);
  auto eg = egads_load(egads_path);
  egads_reclassify(&mesh, eg);
  run_case(&mesh, eg, vtk_path);
  egads_free(eg);
  meshb::write(&mesh, out_path, 2);
  return 0;
}
