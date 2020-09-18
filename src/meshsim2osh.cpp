#include <iostream>

#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

#include "Omega_h_array_ops.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh-in");
  cmdline.add_arg<std::string>("model-in");
  cmdline.add_arg<std::string>("mesh-out");
  cmdline.add_arg<std::string>("model-out");
  if (!cmdline.parse_final(comm, &argc, argv)) return -1;
  auto mesh_in = cmdline.get<std::string>("mesh-in");
  auto model_in = cmdline.get<std::string>("model-in");
  auto mesh_out = cmdline.get<std::string>("mesh-out");
  auto model_out = cmdline.get<std::string>("model-out");
  auto mesh = Omega_h::meshsim::read(mesh_in, model_in, comm);
  Omega_h::binary::write(mesh_out, &mesh);
  std::cout << "wrote mesh " << mesh_out << "\n";

  Omega_h::binary::write_model(model_out, &mesh);
  auto old_model_ents_v = mesh.ask_model_ents(0);
  auto old_model_matches_v = mesh.ask_model_matches(0);
  auto old_model_ents_e = mesh.ask_model_ents(1);
  auto old_model_matches_e = mesh.ask_model_matches(1);
  auto old_model_ents_f = mesh.ask_model_ents(2);
  auto old_model_matches_f = mesh.ask_model_matches(2);
  Omega_h::binary::read_model(model_out, &mesh);
  auto new_model_ents_v = mesh.ask_model_ents(0);
  auto new_model_matches_v = mesh.ask_model_matches(0);
  auto new_model_ents_e = mesh.ask_model_ents(1);
  auto new_model_matches_e = mesh.ask_model_matches(1);
  auto new_model_ents_f = mesh.ask_model_ents(2);
  auto new_model_matches_f = mesh.ask_model_matches(2);
  OMEGA_H_CHECK(old_model_ents_v == new_model_ents_v);
  OMEGA_H_CHECK(old_model_matches_v == new_model_matches_v);
  OMEGA_H_CHECK(old_model_ents_e == new_model_ents_e);
  OMEGA_H_CHECK(old_model_matches_e == new_model_matches_e);
  OMEGA_H_CHECK(old_model_ents_f == new_model_ents_f);
  OMEGA_H_CHECK(old_model_matches_f == new_model_matches_f);
  return 0;
}
