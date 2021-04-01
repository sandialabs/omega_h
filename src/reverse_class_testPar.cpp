#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_array_ops.hpp"
using namespace Omega_h;

void call_print(LOs a) {
  fprintf(stderr,"\n");
  auto a_w = Write<LO> (a.size());
  auto r2w = OMEGA_H_LAMBDA(LO i) {
    a_w[i] = a[i];
  };
  parallel_for(a.size(), r2w);
  auto a_host = HostWrite<LO>(a_w);
  for (int i=0; i<a_host.size(); ++i) {
    fprintf(stderr," %d,", a_host[i]);
  };
  fprintf(stderr,"\n");
  fprintf(stderr,"\n");
  return;
}

void test_3d_p(Library *lib) {

  auto mesh = Mesh(lib);
  binary::read ("./../../omega_h/meshes/unitbox_cutTriCube_1k_4p.osh",
                lib->world(), &mesh);

  // test reverse class APIs
  OMEGA_H_CHECK (!mesh.has_revClass(3));
  auto reg_rc = mesh.ask_revClass(3);
  auto reg_rc_get = mesh.get_revClass(3);
  OMEGA_H_CHECK (mesh.has_revClass(3));
  OMEGA_H_CHECK (reg_rc.ab2b == reg_rc_get.ab2b);
  OMEGA_H_CHECK (reg_rc.a2ab == reg_rc_get.a2ab);

  OMEGA_H_CHECK (!mesh.has_revClass(2));
  auto face_rc = mesh.ask_revClass(2);
  auto face_rc_get = mesh.get_revClass(2);
  OMEGA_H_CHECK (mesh.has_revClass(2));
  OMEGA_H_CHECK (face_rc.ab2b == face_rc_get.ab2b);
  OMEGA_H_CHECK (face_rc.a2ab == face_rc_get.a2ab);

  OMEGA_H_CHECK (!mesh.has_revClass(1));
  auto edge_rc = mesh.ask_revClass(1);
  auto edge_rc_get = mesh.get_revClass(1);
  OMEGA_H_CHECK (mesh.has_revClass(1));
  OMEGA_H_CHECK (edge_rc.ab2b == edge_rc_get.ab2b);
  OMEGA_H_CHECK (edge_rc.a2ab == edge_rc_get.a2ab);

  OMEGA_H_CHECK (!mesh.has_revClass(0));
  auto vert_rc = mesh.ask_revClass(0);
  auto vert_rc_get = mesh.get_revClass(0);
  OMEGA_H_CHECK (mesh.has_revClass(0));
  OMEGA_H_CHECK (vert_rc.ab2b == vert_rc_get.ab2b);
  OMEGA_H_CHECK (vert_rc.a2ab == vert_rc_get.a2ab);

  return;
}

int main(int argc, char** argv) {

  OMEGA_H_CHECK(argc != -1);
  OMEGA_H_CHECK(std::string(argv[0]) != "");

  auto lib = Library();

  test_3d_p(&lib);
  // using mfem adapt tests, it was confirmed that rc info is
  // destroyed during adapt

  return 0;
}
