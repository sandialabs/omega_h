#include "Omega_h_library.hpp"
#include "Omega_h_teuchos.hpp"
#include "Omega_h_cmdline.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  auto cmdline = Omega_h::CmdLine();
  cmdline.add_arg<std::string>("input.{xml,yaml}", "input config file");
  cmdline.add_arg<std::string>("output.{xml,yaml}", "ouput config file");
  if (!cmdline.parse_all_or_help(comm, &argc, argv)) {
    return 2;
  }
  auto inpath = cmdline.get<std::string>("input.{xml,yaml}");
  auto outpath = cmdline.get<std::string>("output.{xml,yaml}");
  Teuchos::ParameterList params;
  auto comm_teuchos = Omega_h::make_teuchos_comm(comm);
  Omega_h::update_parameters_from_file(inpath, &params, *comm_teuchos);
  Omega_h::write_parameters(outpath, params);
}
