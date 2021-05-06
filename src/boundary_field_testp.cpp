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
static void set_target_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto target_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.1;
    h[dim - 1] = 0.001 + 0.198 * std::abs(z - 0.5);
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(target_metrics_w, v, m);
  };
  parallel_for(mesh->nverts(), f);
  mesh->set_tag(VERT, "target_metric", Reals(target_metrics_w));
}

template <Int dim>
void run_case(Mesh* mesh, char const* vtk_path) {
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
  printf("run case 1\n");
  if (vtk_path) {
    writer = vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }

  printf("run case 2\n");
  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  //opts.xfer_opts.type_map["field1_boundary"] = OMEGA_H_LINEAR_INTERP;
  //opts.xfer_opts.type_map["field2_boundary"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["field3_boundary"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["field4_boundary"] = OMEGA_H_POINTWISE;
  //TODO: change this so that the type map should only be "field1"
  Now t0 = now();
  while (approach_metric(mesh, opts)) {
    printf("adapt 1 \n");
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
    if (vtk_path) writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";

}

void test_3d(Library *lib) {

  auto mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/box_3d_2p.osh",
                lib->world(), &mesh);

  auto vtx_rc = mesh.ask_revClass(0);
  auto vert_boundary_ids = (mesh.ask_revClass(0)).ab2b;
  auto nbvert = vert_boundary_ids.size();
  auto face_boundary_ids = (mesh.ask_revClass(2)).ab2b;
  auto nbface = face_boundary_ids.size();
  auto reg_boundary_ids = (mesh.ask_revClass(3)).ab2b;
  auto nbreg = reg_boundary_ids.size();

  mesh.add_boundaryField<Real>(0, "field2", 1);//field1 is present in input mesh
  mesh.add_boundaryField<Real>(2, "field3", 1);
  mesh.add_boundaryField<Real>(3, "field4", 1);
  Write<Real> valsf(nbface, 10);
  Read<Real> valsf_r(valsf);
  Write<Real> valsr(nbreg, 10);
  Read<Real> valsr_r(valsr);
  mesh.set_boundaryField_array<Real>(2, "field3", valsf_r);
  mesh.set_boundaryField_array<Real>(3, "field4", valsr_r);
  const auto rank = lib->world()->rank();
  if (!rank) {
    Write<Real> vals(nbvert, 100);
    Read<Real> vals_r(vals);
    mesh.set_boundaryField_array<Real>(0, "field2", vals_r);
  }
  else {
    Write<Real> vals(nbvert, 50.45632);
    Read<Real> vals_r(vals);
    mesh.set_boundaryField_array<Real>(0, "field2", vals_r);
  }

  vtk::write_parallel ("./../../omega_h/meshes/box_3d_2p.vtk",
                   &mesh);

  auto new_mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/box_3d_2p.osh",
                lib->world(), &new_mesh);
  auto new_bField = new_mesh.get_boundaryField_array<Real>(0, "field1");
  auto vals_r = mesh.get_boundaryField_array<Real>(0, "field1");
  OMEGA_H_CHECK(new_bField == vals_r);
  auto nverts = mesh.nverts();
  OMEGA_H_CHECK(new_bField.size() <= nverts);

  mesh.sync_boundaryField(0, "field2");
  vtk::write_parallel
    ("./../../omega_h/meshes/box_3d_2p_sync.vtk", &mesh);
  mesh.reduce_boundaryField(0, "field2", OMEGA_H_SUM);
  vtk::write_parallel
    ("./../../omega_h/meshes/box_3d_2p_reduce.vtk", &mesh);

  vtk::write_parallel ("./../../omega_h/meshes/box_3d_foo_2p.vtk",
                   &mesh);
  mesh.remove_boundaryField(0, "field1");
  mesh.remove_boundaryField(0, "field2");
  run_case<3>(&mesh, "./../../omega_h/meshes/adapt/box3d_2p.vtk");

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  test_3d(&lib);

  return 0;
}
