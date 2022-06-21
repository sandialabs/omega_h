#include <iostream>

#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

#include "Omega_h_build.hpp"
int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh-in");
  cmdline.add_arg<std::string>("model-in(geomSim)");
  cmdline.add_arg<std::string>("mesh-out");
  cmdline.add_arg<std::string>("numbering-in");
  if (!cmdline.parse_final(comm, &argc, argv)) return -1;
  auto mesh_in = cmdline.get<std::string>("mesh-in");
  auto model_in = cmdline.get<std::string>("model-in(geomSim)");
  auto mesh_out = cmdline.get<std::string>("mesh-out");
  auto numbering_in = cmdline.get<std::string>("numbering-in");
  auto mesh = Omega_h::meshsim::read(mesh_in, model_in, numbering_in, comm);
  auto family = mesh.family();
  if ((family == OMEGA_H_SIMPLEX) || family == OMEGA_H_HYPERCUBE) {
    Omega_h::binary::write(mesh_out, &mesh);
  }

  return 0;
}
