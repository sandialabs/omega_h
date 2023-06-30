#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_vtk.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.meshb");
  cmdline.add_arg<std::string>("mach_in.solb");
  cmdline.add_arg<std::string>("metric_in.solb");
#ifdef OMEGA_H_USE_EGADS
  cmdline.add_arg<std::string>("geom_in.egads");
#endif
  cmdline.add_arg<std::string>("out_prefix");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto mesh_path = cmdline.get<std::string>("mesh_in.meshb");
  auto mach_path = cmdline.get<std::string>("mach_in.solb");
  auto metric_path = cmdline.get<std::string>("metric_in.solb");
#ifdef OMEGA_H_USE_EGADS
  auto geom_path = cmdline.get<std::string>("geom_in.egads");
#endif
  auto out_prefix = cmdline.get<std::string>("out_prefix");
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, mesh_path.c_str());
  Omega_h::meshb::read_sol(&mesh, mach_path.c_str(), "mach");
  Omega_h::meshb::read_sol(&mesh, metric_path.c_str(), "original_metric");
#ifdef OMEGA_H_USE_EGADS
  auto geom = Omega_h::egads_load(geom_path);
#endif
  auto original_metrics =
      mesh.get_array<Omega_h::Real>(Omega_h::VERT, "original_metric");
  fprintf(stderr, "nverts %d nelms %d\n", mesh.nverts(), mesh.nelems());
  fprintf(stderr, "original_metrics size %d\n", original_metrics.size());
  auto graded_metrics =
      Omega_h::limit_metric_gradation(&mesh, original_metrics, 1.0);
  mesh.add_tag(Omega_h::VERT, "target_metric", Omega_h::symm_ncomps(mesh.dim()),
      graded_metrics);
  Omega_h::binary::write("deltaWing.osh", &mesh);
  Omega_h::vtk::write_parallel("deltaWing.vtk", &mesh);
  Omega_h::add_implied_metric_tag(&mesh);
  mesh.ask_qualities();
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.xfer_opts.type_map["mach"] = OMEGA_H_LINEAR_INTERP;
  opts.min_quality_allowed = 0.1;
#ifdef OMEGA_H_USE_EGADS
  opts.egads_model = geom;
#endif

  auto vtx_rc = mesh.ask_revClass(0);
  auto vert_rc_ids = (mesh.ask_revClass(0)).ab2b;
  auto nbvert = vert_rc_ids.size();
  OMEGA_H_CHECK (nbvert < mesh.nverts());
  auto edge_rc_ids = (mesh.ask_revClass(1)).ab2b;
  auto nbedge = edge_rc_ids.size();
  auto face_rc = mesh.ask_revClass(2);
  auto face_a2abSize = face_rc.a2ab.size();
  OMEGA_H_CHECK(face_a2abSize);
  auto face_rc_ids = (mesh.ask_revClass(2)).ab2b;
  auto nbface = face_rc_ids.size();
  auto reg_rc_ids = (mesh.ask_revClass(3)).ab2b;
  auto nbreg = reg_rc_ids.size();

  mesh.add_rcField<Omega_h::LO>(0, "field", 1);
  Omega_h::Write<Omega_h::LO> valsr_v(nbvert, 100);
  mesh.set_rcField_array(0, "field", Omega_h::Read<Omega_h::LO>(valsr_v));
  Omega_h::Write<Omega_h::Real> vals_real(nbvert, 50.2523232);
  mesh.add_rcField<Omega_h::Real>(0, "field1", 1, Omega_h::Read<Omega_h::Real>(vals_real));
  mesh.add_rcField<Omega_h::LO>(1, "field", 1);
  Omega_h::Write<Omega_h::LO> vals(nbedge, 100);
  Omega_h::Read<Omega_h::LO> vals_r(vals);
  mesh.set_rcField_array(1, "field", vals_r);
  mesh.add_rcField<Omega_h::LO>(2, "field", 1);
  Omega_h::Write<Omega_h::LO> valsf(nbface, 12);
  mesh.set_rcField_array(2, "field", Omega_h::Read<Omega_h::LO>(valsf));
  mesh.add_rcField<Omega_h::LO>(3, "field", 1);
  Omega_h::Write<Omega_h::LO> valsr(nbreg, 100);
  Omega_h::Read<Omega_h::LO> valsr_r(valsr);
  mesh.set_rcField_array(3, "field", valsr_r);

  add_rcField_transferMap(&opts, "field", OMEGA_H_INHERIT);
  add_rcField_transferMap(&opts, "field1", OMEGA_H_LINEAR_INTERP);

  opts.verbosity = Omega_h::EXTRA_STATS;
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
  }
#ifdef OMEGA_H_USE_EGADS
  Omega_h::egads_free(geom);
#endif
  Omega_h::vtk::write_parallel("deltaWing_adapted_omegah.vtk", &mesh);
  return 0;
}
