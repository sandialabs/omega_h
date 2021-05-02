#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_array_ops.hpp"

using namespace Omega_h;

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

  printf("btest ok1\n");
  vtk::write_vtu ("./../../omega_h/meshes/plate_6elem.vtu",
                   &mesh);
  printf("btest ok2\n");
  binary::write ("./../../omega_h/meshes/plate_6elem_bField.osh",
                   &mesh);
  printf("btest ok3\n");

  auto new_mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/plate_6elem_bField.osh",
                lib->world(), &new_mesh);
  printf("btest ok4\n");
  auto new_bField = new_mesh.get_boundaryField_array<Real>(0, "field1"); 
  printf("btest ok5\n");
  OMEGA_H_CHECK(new_bField == vals_r);
  OMEGA_H_CHECK(new_bField.size() < nverts);

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  test_2d(&lib);

  return 0;
}
