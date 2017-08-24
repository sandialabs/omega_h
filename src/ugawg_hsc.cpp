#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_vtk.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.meshb");
  cmdline.add_arg<std::string>("mach_in.solb");
  cmdline.add_arg<std::string>("metric_in.solb");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto mesh_path = cmdline.get<std::string>("mesh_in.meshb");
  auto mach_path = cmdline.get<std::string>("mach_in.solb");
  auto metric_path = cmdline.get<std::string>("metric_in.solb");
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, mesh_path.c_str());
  Omega_h::meshb::read_sol(&mesh, mach_path.c_str(), "mach");
  Omega_h::meshb::read_sol(&mesh, metric_path.c_str(), "original_metric");
  auto original_metrics = mesh.get_array<Omega_h::Real>(Omega_h::VERT, "original_metric");
  auto graded_metrics = Omega_h::limit_metric_gradation(&mesh,
      original_metrics, 1.0);
  mesh.add_tag(Omega_h::VERT, "target_metric",
      Omega_h::symm_ncomps(mesh.dim()), graded_metrics);
  Omega_h::add_implied_metric_tag(&mesh);
  mesh.ask_qualities();
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.xfer_opts.type_map["mach"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["original_metric"] = OMEGA_H_METRIC;
  Omega_h::vtk::Writer writer("ugawg_hsc_vtk", &mesh);
  writer.write();
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
    writer.write();
  }
  return 0;
}

