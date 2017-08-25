#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_vtk.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif

#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.meshb");
  cmdline.add_arg<std::string>("mach_in.solb");
  cmdline.add_arg<std::string>("metric_in.solb");
  cmdline.add_arg<std::string>("geom_in.egads");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto mesh_path = cmdline.get<std::string>("mesh_in.meshb");
  auto mach_path = cmdline.get<std::string>("mach_in.solb");
  auto metric_path = cmdline.get<std::string>("metric_in.solb");
#ifdef OMEGA_H_USE_EGADS
  auto geom_path = cmdline.get<std::string>("geom_in.egads");
#endif
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, mesh_path.c_str());
  Omega_h::vtk::write_vtu("classified_faces.vtu", &mesh, 1);
  Omega_h::vtk::write_vtu("classified_edges.vtu", &mesh, 1);
  Omega_h::meshb::read_sol(&mesh, mach_path.c_str(), "mach");
  Omega_h::meshb::read_sol(&mesh, metric_path.c_str(), "original_metric");
#ifdef OMEGA_H_USE_EGADS
  auto geom = Omega_h::egads_load(geom_path);
  Omega_h::egads_reclassify(&mesh, geom);
  Omega_h::vtk::write_vtu("reclassified_faces.vtu", &mesh, 1);
  Omega_h::vtk::write_vtu("reclassified_edges.vtu", &mesh, 1);
#endif
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
#ifdef OMEGA_H_USE_EGADS
  opts.egads_model = geom;
#endif
  Omega_h::vtk::Writer writer("ugawg_hsc_vtk", &mesh);
  writer.write();
  std::cout << "wrote!\n";
  lib.world()->barrier();
  while (Omega_h::approach_metric(&mesh, opts)) {
    std::cout << "before adapt!\n";
    Omega_h::adapt(&mesh, opts);
    std::cout << "after adapt!\n";
    writer.write();
    std::cout << "wrote again!\n";
    lib.world()->barrier();
  }
#ifdef OMEGA_H_USE_EGADS
  Omega_h::egads_free(geom);
#endif
  return 0;
}
