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
  if (!cmdline.parse_final(comm, &argc, argv)) return -1;
  auto mesh_in = cmdline.get<std::string>("mesh-in");
  auto model_in = cmdline.get<std::string>("model-in(geomSim)");
  auto mesh_out = cmdline.get<std::string>("mesh-out");
  printf(" ok1 \n");
  auto mesh = Omega_h::meshsim::read(mesh_in, model_in, comm);
  printf(" ok2 \n");
  //Omega_h::binary::write(mesh_out, &mesh);
  //std::cout << "wrote mesh " << mesh_out << "\n";

  //testing for codes generation
  //auto mesh_build = Omega_h::build_box(lib.world(), OMEGA_H_HYPERCUBE, 1.0, 1.0, 0.0, 2, 2, 0);


  return 0;
}