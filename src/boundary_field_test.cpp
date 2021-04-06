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
  mesh.add_tag<Real>(0, "field", 3);
  Write<Real> vals(nvert*3, 50, "fvals");
  Read<Real> vals_r(vals);
  mesh.set_tag<Real>(0, "field", vals_r);
  mesh.change_tagToBoundary<Real>(0, 3, "field");
  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);

  test_2d(&lib);

  return 0;
}
