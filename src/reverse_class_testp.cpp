#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_array_ops.hpp"
using namespace Omega_h;

void test_cutBox(Library *lib, const std::string &mesh_file) {

  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);

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

  auto rc_r2f = mesh.ask_revClass_downAdj (3, 2);
  auto rc_r2e = mesh.ask_revClass_downAdj (3, 1);
  auto rc_r2v = mesh.ask_revClass_downAdj (3, 0);
  auto rc_face2vert = mesh.ask_revClass_downAdj (2, 0);
  auto rc_face2e = mesh.ask_revClass_downAdj (2, 1);
  auto rc_e2v = mesh.ask_revClass_downAdj (1, 0);

  return;
}

void test_box(Library *lib, const std::string &mesh_file) {
 
  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);
  auto rank = lib->world()->rank();
 
  if (!rank) {
    auto vtx2_6_rc = mesh.ask_revClass(0, LOs({2, 6}));
    OMEGA_H_CHECK (vtx2_6_rc.ab2b == LOs({17}));
    OMEGA_H_CHECK (vtx2_6_rc.a2ab == LOs({0, 1, 1}));
  }

  if (rank) {
    auto vtx2_6_rc = mesh.ask_revClass(0, LOs({2, 6}));
    OMEGA_H_CHECK (vtx2_6_rc.ab2b == LOs({7}));
    OMEGA_H_CHECK (vtx2_6_rc.a2ab == LOs({0, 0, 1}));

    auto vtx6_rc = mesh.ask_revClass(0, LOs({6}));
    OMEGA_H_CHECK (vtx6_rc.ab2b == LOs({7}));
    OMEGA_H_CHECK (vtx6_rc.a2ab == LOs({0, 1}));
  }

  return;
}

int main(int argc, char** argv) {

  auto lib = Library(&argc, &argv);

  if (argc != 3) {
    Omega_h_fail("a.out <3d_in_boxMesh> <3d_in_cutMesh>\n");
  }
  char const* path_box = nullptr;
  char const* path_cutBox = nullptr;
  path_box = argv[1];
  path_cutBox = argv[2];

  test_cutBox(&lib, path_cutBox);
  test_box(&lib, path_box);

  return 0;
}
