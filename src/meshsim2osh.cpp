#include <iostream>

#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

#include "Omega_h_array_ops.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();

  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh-in");
  cmdline.add_arg<std::string>("model-in");
  cmdline.add_arg<std::string>("mesh-out");
  cmdline.add_arg<std::string>("model-out");

  cmdline.add_arg<int>("nparts_out");

  auto nparts_total = world->size();
  auto nparts_in = 1;
  //todo: add checks on nparts out and total

  if (!cmdline.parse_final(world, &argc, argv)) return -1;
  auto mesh_in = cmdline.get<std::string>("mesh-in");
  auto model_in = cmdline.get<std::string>("model-in");
  auto mesh_out = cmdline.get<std::string>("mesh-out");
  auto model_out = cmdline.get<std::string>("model-out");

  auto nparts_out = cmdline.get<int>("nparts_out");
  auto is_in = (world->rank() < nparts_in);
  auto comm_in = world->split(int(is_in), 0);
  auto is_out = (world->rank() < nparts_out);
  auto comm_out = world->split(int(is_out), 0);

  if (is_in) {
  //if (!world->rank()) {
  printf("nparts1out=%d\n", nparts_out);
    auto mesh = Omega_h::meshsim::read(mesh_in, model_in, comm_in);
  printf("nparts2out=%d\n", nparts_out);
    Omega_h::binary::write(mesh_out, &mesh);
  printf("nparts3out=%d\n", nparts_out);
    Omega_h::binary::write_model(model_out, &mesh);
  printf("nparts4out=%d\n", nparts_out);
  }

  return 0;
}
