#include "Omega_h.hpp"
#include "Omega_h_egads.hpp"
#include "Omega_h_math.hpp"
#include "loop.hpp"
#include "metric.hpp"
#include "timer.hpp"

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 3;

static Reals get_first_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_dofs(dim));
  constexpr Real h0 = 0.001;
  auto f = LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto z = p[2];
    auto h = vector_3(0.1, 0.1, h0 + 2 * (0.1 - h0) * fabs(z - 0.5));
    auto m = diagonal(metric_eigenvalues(h));
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_second_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_dofs(dim));
  constexpr Real h0 = 0.001;
  // constexpr Real h0 = 0.004;
  auto f = LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = atan2(y, x);
    auto h = vector_3(h0 + 2 * (0.1 - h0) * fabs(radius - 0.5), 0.1, 0.1);
    auto rotation = rotate(t, vector_3(0, 0, 1));
    auto m = compose_metric(rotation, h);
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static void set_target_metric(Mesh* mesh, int which_metric) {
  auto target_metrics =
      (which_metric == 1) ? get_first_metric(mesh) : get_second_metric(mesh);
  mesh->set_tag(VERT, "target_metric", target_metrics);
}

static void run_case(
    Mesh* mesh, Egads* eg, int which_metric, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = find_implied_metric(mesh);
  mesh->add_tag(VERT, "metric", symm_dofs(dim), OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT, implied_metrics);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->add_tag<Real>(
      VERT, "target_metric", symm_dofs(dim), OMEGA_H_METRIC, OMEGA_H_DO_OUTPUT);
  set_target_metric(mesh, which_metric);
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
  opts.egads_model = eg;
  if (which_metric == 2) {
    opts.should_move_for_quality = true;
    opts.max_length_allowed = opts.max_length_desired * 2.0;
    opts.min_quality_desired = 0.25;
    opts.min_quality_allowed = 0.20;
  } else {
    opts.max_length_allowed = opts.max_length_desired * 2.0;
  }
  Now t0 = now();
  while (approach_size_field(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) {
      set_target_metric(mesh, which_metric);
    }
    if (vtk_path) writer.write();
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
  int which_metric = 1;
  bool should_help = false;
  for (int i = 1; i < argc; ++i) {
    if (std::string("-a") == argv[i]) {
      if (i == argc - 1) {
        std::cout << "-a takes an argument\n";
        should_help = true;
      } else {
        vtk_path = argv[++i];
      }
    } else if (std::string("-m") == argv[i]) {
      if (i == argc - 1) {
        std::cout << "-m takes an argument\n";
        should_help = true;
      } else {
        which_metric = atoi(argv[++i]);
        if (!(which_metric == 1 || which_metric == 2)) {
          std::cout << "-m takes only values 1 or 2\n";
          should_help = true;
        }
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
    std::cout
        << "options: -m (1|2)                   1 for shock, 2 for BL\n";
    return -1;
  }
  Mesh mesh(&lib);
  meshb::read(&mesh, in_path);
  auto eg = egads_load(egads_path);
  egads_reclassify(&mesh, eg);
  run_case(&mesh, eg, which_metric, vtk_path);
  egads_free(eg);
  meshb::write(&mesh, out_path, 2);
  return 0;
}
