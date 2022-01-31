#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_array_ops.hpp"
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_build.hpp>

using namespace Omega_h;

template <Int dim>
void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  auto calc_metrics = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.5;
    h[dim - 1] = 2*(0.001 + 0.198 * std::abs(z - 0.5));
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), std::move(calc_metrics), "calc_metrics");
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
}

template <Int dim>
void run_adapt(Mesh* mesh, char const* vtk_path) {
  auto world = mesh->comm();
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto implied_metrics = get_implied_metrics(mesh);
  mesh->add_tag(0, "metric", symm_ncomps(dim), implied_metrics);
  mesh->add_tag<Real>(VERT, "target_metric", symm_ncomps(dim));
  set_target_metric<dim>(mesh);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
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
}

void test_3d(Library *lib, const std::string &mesh_file, const char* vtk_file,
             const char *sync_file, const char *reduce_file,
             const char *adapt_file) {

  // TODO: change the adapt b.field transfer test to work for FACE
  // without having to add pseudo rcFields to all other dims
  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);

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

  mesh.add_rcField<LO>(0, "field", 1);
  Write<LO> valsr_v(nbvert, 100);
  mesh.set_rcField_array(0, "field", Read<LO>(valsr_v));

  mesh.add_rcField<LO>(1, "field", 1);
  const auto rank = lib->world()->rank();
  if ((!rank)) {
    Write<LO> vals(nbedge, 100);
    Read<LO> vals_r(vals);
    mesh.set_rcField_array(1, "field", vals_r);
  }
  else {
    Write<LO> vals(nbedge, 50.45632);
    Read<LO> vals_r(vals);
    mesh.set_rcField_array(1, "field", vals_r);
  }

  mesh.add_rcField<LO>(2, "field", 1);
  Write<LO> valsf(nbface, 12);
  auto ab2b = face_rc.ab2b;
  auto a2ab = face_rc.a2ab;
  auto assign_field_vals = OMEGA_H_LAMBDA(LO gf) {
    if ((gf == 16) || (gf == 22)) {
      auto start = a2ab[gf];
      auto end = a2ab[gf+1];
      for (int index = start; index < end; ++index) {
        valsf[index] = gf;
      }
    }
  };
  parallel_for(face_a2abSize-1, std::move(assign_field_vals),
               "assign_field_vals");
  mesh.set_rcField_array(2, "field", Read<LO>(valsf));
 
  mesh.add_rcField<LO>(3, "field", 1);
  Write<LO> valsr(nbreg, 100);
  Read<LO> valsr_r(valsr);
  mesh.set_rcField_array(3, "field", valsr_r);

  vtk::write_parallel (vtk_file, &mesh);

  mesh.sync_rcField(1, "field");
  vtk::write_parallel (sync_file, &mesh);
  mesh.reduce_rcField(1, "field", OMEGA_H_SUM);
  vtk::write_parallel (reduce_file, &mesh);

  run_adapt<3>(&mesh, adapt_file);

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);
  auto world = lib.world();

  if (argc != 6) {
    Omega_h_fail("a.out <in_mesh> <out_vtk> <out_sync_vtk> <out_reduce_vtk> <out_adapt>\n");
  };
  char const* path_in = nullptr;
  char const* path_vtk = nullptr;
  char const* path_sync = nullptr;
  char const* path_reduce = nullptr;
  char const* path_adapt = nullptr;
  path_in = argv[1];
  path_vtk = argv[2];
  path_sync = argv[3];
  path_reduce = argv[4];
  path_adapt = argv[5];

  world->barrier();
  test_3d(&lib, path_in, path_vtk, path_sync, path_reduce, path_adapt);
  world->barrier();

}
