#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_array_ops.hpp"
#include <Omega_h_metric.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_timer.hpp>
#include <Omega_h_for.hpp>

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
  if (vtk_path) {
    writer = vtk::FullWriter(vtk_path, mesh);
    writer.write();
  }

  auto opts = AdaptOpts(mesh);
  opts.verbosity = EXTRA_STATS;
  opts.length_histogram_max = 2.0;
  opts.max_length_allowed = opts.max_length_desired * 2.0;
  add_boundaryField_transferMap(&opts, "field1", OMEGA_H_LINEAR_INTERP);
  Now t0 = now();
  while (approach_metric(mesh, opts)) {
    adapt(mesh, opts);
    if (mesh->has_tag(VERT, "target_metric")) set_target_metric<dim>(mesh);
    if (vtk_path) writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";

}
void test_2d(Library *lib) {

  auto mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/plate_6elem.osh",
                lib->world(), &mesh);

  auto nverts = mesh.nverts();
  auto boundary_ids = (mesh.ask_revClass(VERT)).ab2b;
  auto nbvert = boundary_ids.size();
  Write<Real> vals(nbvert, 50);
  Read<Real> vals_r(vals);

  mesh.add_boundaryField<Real>(0, "field1", 1, vals_r);

  Write<Real> vals2(mesh.nfaces()*2, 50);
  Read<Real> vals_r2(vals2);

  Write<Real> vals3(mesh.nverts()*3, 500);
  Read<Real> vals_r3(vals3);

  vtk::write_vtu ("./../../omega_h/meshes/plate_6elem.vtu",
                   &mesh);
  binary::write ("./../../omega_h/meshes/plate_6elem_bField.osh",
                   &mesh);

  auto new_mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/plate_6elem_bField.osh",
                lib->world(), &new_mesh);
  auto new_bField = new_mesh.get_boundaryField_array<Real>(0, "field1"); 
  OMEGA_H_CHECK(new_bField == vals_r);
  OMEGA_H_CHECK(new_bField.size() < nverts);

  return;
}

void test_3d(Library *lib) {
  auto mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/box_3d.osh",
                lib->world(), &mesh);
  OMEGA_H_CHECK(mesh.has_boundaryField(0, "field1"));
  run_case<3>(&mesh, "./../../omega_h/meshes/adapt/box3d.vtk");

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  test_2d(&lib);
  test_3d(&lib);

  return 0;
}
