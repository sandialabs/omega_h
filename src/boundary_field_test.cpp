#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
using namespace Omega_h;

void test_2d(Library *lib) {

  auto mesh = Mesh(lib);
  binary::read ("/lore/joshia5/Meshes/oh-mfem/plate_6elem.osh",
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

  OMEGA_H_CHECK(argc != -1);
  OMEGA_H_CHECK(std::string(argv[0]) != "");

  auto lib = Library();

  test_2d(&lib);

  return 0;
}
