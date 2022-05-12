#include <Omega_h_adapt.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>
#include <iostream>
#include <string>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

using namespace Omega_h;

template <Int dim>
void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.5;
    h[dim - 1] = 0.001 + 0.198 * std::abs(z - 0.5);
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
  return;
}

template <Int dim>
void run_adapt(Mesh* mesh, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(0, "metric", symm_ncomps(dim), implied_metrics);
  mesh->add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  set_target_metric<dim>(mesh);

  OMEGA_H_CHECK(mesh->has_tag(0, "metric"));
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  OMEGA_H_CHECK(mesh->has_tag(0, "metric"));

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
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  add_rcField_transferMap(&opts, "field", OMEGA_H_INHERIT);
  Now t0 = now();
  while (approach_metric(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
    if (vtk_path) writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";

  return;
}

void test_2d(Library* lib, const std::string& mesh_file, const char* vtu_file,
    const char* rcField_file) {
  auto mesh = Mesh(lib);
  binary::read(mesh_file, lib->world(), &mesh);

  auto nverts = mesh.nverts();
  auto rc_ids = (mesh.ask_revClass(VERT)).ab2b;
  auto nbvert = rc_ids.size();
  Write<Real> vals(nbvert, 50);
  Read<Real> vals_r(vals);
  mesh.add_rcField<Real>(0, "field1", 1, vals_r);

  auto face_1_2_rc = mesh.ask_revClass(2, LOs({1, 2}));
  auto face_1_2_rc_ids = face_1_2_rc.ab2b;
  Write<Real> vals_f_1_2(face_1_2_rc_ids.size(), 50);
  mesh.add_rcField<Real>(LOs({1, 2}), 2, "face_1_2", 1);
  mesh.set_rcField_array(2, "face_1_2", Reals(vals_f_1_2));

  Write<Real> vals2(mesh.nfaces() * 2, 50);
  Read<Real> vals_r2(vals2);

  Write<Real> vals3(mesh.nverts() * 3, 500);
  Read<Real> vals_r3(vals3);

  vtk::write_vtu(vtu_file, &mesh);
  binary::write(rcField_file, &mesh);

  auto new_mesh = Mesh(lib);
  binary::read(rcField_file, lib->world(), &new_mesh);

  return;
}

void test_3d(Library* lib, const std::string& mesh_file, const char* out_file) {
  auto mesh = Mesh(lib);
  binary::read(mesh_file, lib->world(), &mesh);

  auto face_22_rc = mesh.ask_revClass(2, LOs({22}));
  auto face_22_rc_ids = face_22_rc.ab2b;
  Write<Real> vals_f_22(face_22_rc_ids.size(), 50);
  mesh.add_rcField<Real>(LOs({22}), 2, "face_22", 1);
  mesh.set_rcField_array(2, "face_22", Reals(vals_f_22));
  auto array0 = mesh.get_rcField_array<Real>(2, "face_22");
  OMEGA_H_CHECK(array0.exists());
  binary::write("box3d_withFace22Field.osh", &mesh);
  auto array = mesh.get_rcField_array<Real>(2, "face_22");
  OMEGA_H_CHECK(array.exists());
  auto new_mesh = Mesh(lib);
  binary::read("box3d_withFace22Field.osh", lib->world(), &new_mesh);
  auto vals_f_22_read = new_mesh.get_rcField_array<Real>(2, "face_22");
  OMEGA_H_CHECK(vals_f_22_read == Reals(vals_f_22));

  auto vert_rc_ids = (mesh.ask_revClass(0)).ab2b;
  auto nbvert = vert_rc_ids.size();
  auto edge_rc_ids = (mesh.ask_revClass(1)).ab2b;
  auto nbedge = edge_rc_ids.size();
  auto face_rc_ids = (mesh.ask_revClass(2)).ab2b;
  auto nbface = face_rc_ids.size();
  auto reg_rc_ids = (mesh.ask_revClass(3)).ab2b;
  auto nbreg = reg_rc_ids.size();

  Write<Real> vals_v(nbvert, 100);
  Write<Real> vals_e(nbedge, 100);
  Write<Real> vals_f(nbface, 100);
  Write<Real> vals_r(nbreg, 100);
  mesh.add_rcField<Real>(0, "field", 1, Reals(vals_v));
  mesh.add_rcField<Real>(1, "field", 1, Reals(vals_e));
  mesh.add_rcField<Real>(2, "field", 1, Reals(vals_f));
  mesh.add_rcField<Real>(3, "field", 1, Reals(vals_r));

  run_adapt<3>(&mesh, out_file);

  return;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();

  if (argc != 6) {
    Omega_h_fail(
        "a.out <2d_in_mesh> <3d_in_mesh> <3d_out_mesh> <2d_out_vtu> "
        "<2d_out_osh>\n");
  };
  char const* path_2d = nullptr;
  char const* path_3d = nullptr;
  char const* path_3d_out = nullptr;
  char const* path_2d_vtu = nullptr;
  char const* path_2d_rcField = nullptr;
  path_2d = argv[1];
  path_3d = argv[2];
  path_3d_out = argv[3];
  path_2d_vtu = argv[4];
  path_2d_rcField = argv[5];

  test_2d(&lib, path_2d, path_2d_vtu, path_2d_rcField);
  world->barrier();
  test_3d(&lib, path_3d, path_3d_out);
  world->barrier();
}
