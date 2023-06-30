#include <iostream>
#include <string>

#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_array_ops.hpp"

using namespace Omega_h;

void test_2d(Library *lib, const std::string &mesh_file) {

  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);

  OMEGA_H_CHECK (!mesh.has_revClass(2));
  auto face_rc = mesh.ask_revClass(2);
  auto face_1_2_rc = mesh.ask_revClass(2, LOs({1, 2}));
  auto face_rc_get = mesh.get_revClass(2);
  OMEGA_H_CHECK (mesh.has_revClass(2));
  OMEGA_H_CHECK (face_rc.ab2b == face_rc_get.ab2b);
  OMEGA_H_CHECK (face_rc.a2ab == face_rc_get.a2ab);
  OMEGA_H_CHECK (face_rc.ab2b == LOs({0, 1, 2, 3, 4, 5}));
  OMEGA_H_CHECK (face_rc.a2ab == LOs({0, 0, 0, 6}));
  OMEGA_H_CHECK (face_1_2_rc.ab2b == LOs({0, 1, 2, 3, 4, 5}));
  OMEGA_H_CHECK (face_1_2_rc.a2ab == LOs({0, 0, 6}));

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
  OMEGA_H_CHECK (vert_rc.ab2b == LOs({3, 2, 1, 0}));
  OMEGA_H_CHECK (vert_rc.a2ab == LOs({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2,
                                      2, 2, 2, 3, 3, 3, 3, 4}));

  auto rc_face2vert = mesh.ask_revClass_downAdj (2, 0);
  auto f_classids = mesh.get_array<ClassId>(2, "class_id");
  OMEGA_H_CHECK (rc_face2vert.ab2b == 
    LOs({4, 7, 3, 7, 4, 5, 4, 2, 5, 7, 5, 6, 6, 0, 7, 5, 1, 6}));
  OMEGA_H_CHECK (rc_face2vert.a2ab == LOs({0, 0, 0, 18}));
  auto rc_face2edge = mesh.ask_revClass_downAdj (2, 1);
  auto rc_e2v = mesh.ask_revClass_downAdj (1, 0);

  return;
}

void test_3d(Library *lib, const std::string &mesh_file) {

  auto mesh = Mesh(lib);
  binary::read (mesh_file, lib->world(), &mesh);

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

  auto rc_face2vert = mesh.ask_revClass_downAdj (2, 0);
  auto rc_face2edge = mesh.ask_revClass_downAdj (2, 1);
  auto rc_edge2vert = mesh.ask_revClass_downAdj (1, 0);

  return;
}

int main(int argc, char** argv) {

  auto lib = Library (&argc, &argv);
  
  if (argc != 3) {
    Omega_h_fail("a.out <2d_in_mesh> <3d_in_mesh>\n");
  }
  char const* path_2d = nullptr;
  char const* path_3d = nullptr;
  path_2d = argv[1];
  path_3d = argv[2];

  test_2d(&lib, path_2d);
  test_3d(&lib, path_3d);

  return 0;
}
