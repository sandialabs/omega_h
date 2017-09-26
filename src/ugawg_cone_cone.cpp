#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_metric.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif

#include <iostream>

static void compute_implied_metric(Omega_h::Mesh* mesh) {
  auto metrics = Omega_h::get_implied_metrics(mesh);
  metrics = Omega_h::limit_metric_gradation(mesh, metrics, 1.0);
  mesh->add_tag(Omega_h::VERT, "metric", Omega_h::symm_ncomps(mesh->dim()), metrics);
}

static void compute_target_metric(Omega_h::Mesh* mesh) {
  auto metric = Omega_h::diagonal(Omega_h::metric_eigenvalues_from_lengths(Omega_h::vector_3(0.1, 0.1, 0.1)));
  auto metrics = Omega_h::repeat_symm(mesh->nverts(), metric);
  mesh->add_tag(Omega_h::VERT, "target_metric", Omega_h::symm_ncomps(mesh->dim()), metrics);
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
  auto& viz_flag = cmdline.add_flag("--viz", "optional VTK progress log");
  viz_flag.add_arg<std::string>("path_vtk");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path_in = cmdline.get<std::string>("mesh_in.osh");
  auto path_out = cmdline.get<std::string>("mesh_out.osh");
  Omega_h::Mesh mesh(&lib);
  std::cout << "reading in " << path_in << '\n';
  Omega_h::binary::read(path_in, lib.world(), &mesh);
  std::cout << "computing metric tags\n";
  compute_implied_metric(&mesh);
  compute_target_metric(&mesh);
  std::cout << "computing minimum quality\n";
  Omega_h::AdaptOpts opts(&mesh);
#ifdef OMEGA_H_USE_EGADS
  auto has_model = cmdline.parsed("--model");
  if (has_model) {
    auto model_path = cmdline.get<std::string>("--model", "model.step");
    std::cout << "reading in " << model_path << '\n';
    auto eg = Omega_h::egads_load(model_path);
    Omega_h::egads_reclassify(&mesh, eg);
    opts.egads_model = eg;
  }
#endif
  auto has_viz = cmdline.parsed("--viz");
  Omega_h::vtk::Writer writer;
  if (has_viz) {
    auto viz_path = cmdline.get<std::string>("--viz", "path_vtk");
    writer = Omega_h::vtk::Writer(viz_path, &mesh);
    writer.write();
  }
  opts.verbosity = Omega_h::EXTRA_STATS;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
    if (has_viz) writer.write();
  }
  std::cout << "writing out " << path_out << '\n';
  mesh.remove_tag(Omega_h::VERT, "metric");
  Omega_h::binary::write(path_out, &mesh);
#ifdef OMEGA_H_USE_EGADS
  if (has_model) {
    Omega_h::egads_free(opts.egads_model);
  }
#endif
}
