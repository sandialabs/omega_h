#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_egads.hpp>

#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.meshb");
  cmdline.add_arg<std::string>("metric_in.solb");
  cmdline.add_arg<std::string>("mesh_out.meshb");
  cmdline.add_arg<std::string>("metric_out.solb");
  auto& model_flag = cmdline.add_flag("--model", "optional EGADS model");
  model_flag.add_arg<std::string>("model.step");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto mesh_path_in = cmdline.get<std::string>("mesh_in.meshb");
  auto mesh_path_out = cmdline.get<std::string>("mesh_out.meshb");
  auto metric_path_in = cmdline.get<std::string>("metric_in.solb");
  auto metric_path_out = cmdline.get<std::string>("metric_out.solb");
  Omega_h::Mesh mesh(&lib);
  std::cout << "reading in " << mesh_path_in << '\n';
  Omega_h::meshb::read(&mesh, mesh_path_in);
  std::cout << "reading in " << metric_path_in << '\n';
  Omega_h::meshb::read_sol(&mesh, metric_path_in, "ugawg_metric");
  std::cout << "writing out loaded.vtu\n";
  Omega_h::vtk::write_vtu("loaded.vtu", &mesh);
  Omega_h::AdaptOpts opts(&mesh);
  opts.xfer_opts.type_map["ugawg_metric"] = OMEGA_H_METRIC;
  if (cmdline.parsed("--model")) {
    auto model_path = cmdline.get<std::string>("--model", "model.step");
    std::cout << "reading in " << model_path << '\n';
    opts.egads_model = Omega_h::egads_load(model_path);
  }
  Omega_h::fix(&mesh, opts, true);
  std::cout << "writing out fixed.vtu\n";
  Omega_h::vtk::write_vtu("fixed.vtu", &mesh);
  std::cout << "limiting gradation on UGAWG metric, setting as target\n";
  auto target_metrics = mesh.get_array<Omega_h::Real>(Omega_h::VERT, "ugawg_metric");
  target_metrics = Omega_h::limit_metric_gradation(&mesh, target_metrics, 1.0);
  mesh.add_tag(Omega_h::VERT, "target_metric", Omega_h::symm_ncomps(mesh.dim()), target_metrics);
  std::cout << "writing out before_adapt.vtu\n";
  Omega_h::vtk::write_vtu("before_adapt.vtu", &mesh);
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
  }
  std::cout << "writing out adapted.vtu\n";
  Omega_h::vtk::write_vtu("adapted.vtu", &mesh);
  std::cout << "writing out " << mesh_path_out << '\n';
  Omega_h::meshb::write(&mesh, mesh_path_out);
  std::cout << "writing out " << metric_path_out << '\n';
  Omega_h::meshb::write_sol(&mesh, metric_path_out, "ugawg_metric");
#ifdef OMEGA_H_USE_EGADS
  if (opts.egads_model != nullptr) {
    Omega_h::egads_free(opts.egads_model);
  }
#endif
}
