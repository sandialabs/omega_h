#include "Omega_h_adapt.hpp"
#include "Omega_h_egads.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_timer.hpp"

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 3;

static Reals get_first_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  constexpr Real h0 = 0.001;
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto z = p[2];
    auto h = vector_3(0.1, 0.1, h0 + 2 * (0.1 - h0) * std::abs(z - 0.5));
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_second_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  constexpr Real h0 = 0.001;
  constexpr Real h_z = 1.0 / 10.0;
  constexpr Real h_t = 1.0 / 10.0;
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = std::atan2(y, x);
    auto h = vector_3(h0 + 2 * (0.1 - h0) * std::abs(radius - 0.5), h_t, h_z);
    auto rotation = rotate(t, vector_3(0, 0, 1));
    auto m = compose_metric(rotation, h);
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_third_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  constexpr Real h0 = 0.001;
  constexpr Real h_z = 1.0 / 10.0;
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = std::atan2(y, x);
    auto d = (0.6 - radius) * 10.0;
    Real h_t = (d < 0.0) ? (1.0 / 10.0)
                         : (d * (1.0 / 40.0) + (1.0 - d) * (1.0 / 10.0));
    auto h = vector_3(h0 + 2 * (0.1 - h0) * std::abs(radius - 0.5), h_t, h_z);
    auto rotation = rotate(t, vector_3(0, 0, 1));
    auto m = compose_metric(rotation, h);
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static void set_target_metric(Mesh* mesh, int which_metric, bool should_limit) {
  auto target_metrics =
      ((which_metric == 1) ? get_first_metric(mesh)
                           : ((which_metric == 2) ? get_second_metric(mesh)
                                                  : get_third_metric(mesh)));
  if (should_limit) {
    target_metrics = limit_metric_gradation(mesh, target_metrics, 1.0);
  }
  mesh->set_tag(VERT, "target_metric", target_metrics);
}

static void run_case(
    Mesh* mesh, Egads* eg, int which_metric, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  add_implied_metric_tag(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  mesh->add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  bool should_limit = (which_metric == 2);
  set_target_metric(mesh, which_metric, should_limit);
  mesh->ask_lengths();
  mesh->ask_qualities();
  vtk::FullWriter writer;
  if (vtk_path) {
    writer = vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }
  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.egads_model = eg;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  Now t0 = now();
  while (approach_metric(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) {
      set_target_metric(mesh, which_metric, should_limit);
    }
    if (vtk_path) writer.write();
  }
  if (which_metric == 3) {
    adapt(mesh, opts);
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
        if (!(which_metric >= 1 && which_metric <= 3)) {
          std::cout << "-m takes only values 1, 2 or 3\n";
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
    std::cout << "options: -a vtk_path      debug output for adaptivity\n";
    std::cout << "options: -m (1|2|3)       1 for shock, 2 for BL, 3 for "
                 "quality BL\n";
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
