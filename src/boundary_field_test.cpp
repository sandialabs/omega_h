#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
using namespace Omega_h;

void test_2d(Library *lib) {

  auto mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/plate_6elem.osh",
                lib->world(), &mesh);

  auto nvert = mesh.nverts();
  mesh.add_tag<Real>(0, "field1", 1);
  Write<Real> vals(nvert, 50);
  Read<Real> vals_r(vals);
  mesh.set_tag<Real>(0, "field1", vals_r);

  mesh.change_tagToBoundary<Real>(0, 1, "field1");
  OMEGA_H_CHECK(!mesh.has_tag(0, "field1"));
  OMEGA_H_CHECK(mesh.has_tag(0, "field1_boundary"));

  mesh.change_tagToMesh<Real>(0, 1, "field1_boundary");
  OMEGA_H_CHECK(mesh.has_tag(0, "field1"));
  OMEGA_H_CHECK(!mesh.has_tag(0, "field1_boundary"));

  mesh.add_tag<Real>(2, "field2", 2);
  Write<Real> vals2(mesh.nfaces()*2, 50);
  Read<Real> vals_r2(vals2);
  mesh.set_tag<Real>(2, "field2", vals_r2);

  mesh.change_tagToBoundary<Real>(2, 2, "field2");
  OMEGA_H_CHECK(!mesh.has_tag(2, "field2"));
  OMEGA_H_CHECK(mesh.has_tag(2, "field2_boundary"));

  mesh.change_tagToMesh<Real>(2, 2, "field2_boundary");
  OMEGA_H_CHECK(mesh.has_tag(2, "field2"));
  OMEGA_H_CHECK(!mesh.has_tag(2, "field2_boundary"));

  mesh.add_tag<Real>(0, "meshfield", 3);
  Write<Real> vals3(mesh.nverts()*3, 500);
  Read<Real> vals_r3(vals3);
  mesh.set_tag<Real>(0, "meshfield", vals_r3);

  vtk::write_vtu ("./../../omega_h/meshes/plate_6elem.vtu",
                   &mesh);
  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  test_2d(&lib);

  return 0;
}
