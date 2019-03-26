#include <iostream>

#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  Omega_h::CmdLine cmdline;
  auto& mesh_in_flag = cmdline.add_flag("--mesh-in", "input mesh file path");
  mesh_in_flag.add_arg<std::string>("mesh.sms");
  auto& model_in_flag = cmdline.add_flag("--model-in", "input model file path");
  model_in_flag.add_arg<std::string>("model.smd");
  if (!cmdline.parse_final(comm, &argc, argv)) {
    return -1;
  }
  if (!cmdline.parsed("--mesh-in")) {
    std::cout << "No input mesh file specified\n";
    cmdline.show_help(argv);
    return -1;
  }
  if (!cmdline.parsed("--model-in")) {
    std::cout << "No input model file specified\n";
    cmdline.show_help(argv);
    return -1;
  }
  Omega_h::filesystem::path mesh_in =
      cmdline.get<std::string>("--mesh-in", "mesh.sms");
  Omega_h::filesystem::path model_in =
      cmdline.get<std::string>("--model-in", "model.smd");
  std::cout << "Loading mesh from " << mesh_in << "\n";
  std::cout << "Loading model from " << model_in << "\n";
  auto mesh = Omega_h::meshsim::read(mesh_in, model_in, comm);
}
