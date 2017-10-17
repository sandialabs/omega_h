#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif

#include <iostream>

static void compute_metric(Omega_h::Mesh* mesh) {
  auto metrics = Omega_h::get_implied_metrics(mesh);
  metrics = Omega_h::limit_metric_gradation(mesh, metrics, 1.0);
  mesh->remove_tag(Omega_h::VERT, "metric");
  mesh->add_tag(
      Omega_h::VERT, "metric", Omega_h::symm_ncomps(mesh->dim()), metrics);
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.osh");
  cmdline.add_arg<std::string>("mesh_out.osh");
#ifdef OMEGA_H_USE_EGADS
  auto& model_flag = cmdline.add_flag("--model", "optional EGADS model");
  model_flag.add_arg<std::string>("model.step");
#endif
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path_in = cmdline.get<std::string>("mesh_in.osh");
  auto path_out = cmdline.get<std::string>("mesh_out.osh");
  Omega_h::Mesh mesh(&lib);
  std::cout << "reading in " << path_in << '\n';
  Omega_h::binary::read(path_in, lib.world(), &mesh);
  std::cout << "computing implied metric tag\n";
  compute_metric(&mesh);
  Omega_h::Real minqual_old = -1.0;
  std::cout << "computing minimum quality\n";
  auto minqual = mesh.min_quality();
  std::cout << "minimum quality " << minqual << '\n';
  Omega_h::AdaptOpts opts(&mesh);
#ifdef OMEGA_H_USE_EGADS
  if (cmdline.parsed("--model")) {
    auto model_path = cmdline.get<std::string>("--model", "model.step");
    std::cout << "reading in " << model_path << '\n';
    auto eg = Omega_h::egads_load(model_path);
    opts.egads_model = eg;
  }
#endif
  while (minqual < opts.min_quality_desired && minqual > minqual_old) {
    std::cout << "current quality " << minqual << '\n';
    minqual_old = minqual;
    opts.min_quality_allowed = minqual;
    opts.verbosity = Omega_h::EXTRA_STATS;
    opts.nsliver_layers = 10;
    Omega_h::adapt(&mesh, opts);
    std::cout << "recomputing implied metric tag\n";
    compute_metric(&mesh);
    std::cout << "recomputing minimum quality\n";
    minqual = mesh.min_quality();
    std::cout << "minimum quality " << minqual << '\n';
  }
  std::cout << "writing out " << path_out << '\n';
  mesh.remove_tag(Omega_h::VERT, "metric");
  Omega_h::binary::write(path_out, &mesh);
#ifdef OMEGA_H_USE_EGADS
  if (opts.egads_model != nullptr) {
    Omega_h::egads_free(opts.egads_model);
  }
#endif
}
