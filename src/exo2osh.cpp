#include <Omega_h.hpp>
#include <Omega_h_cmdline.hpp>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  CmdLine cmdline;
  cmdline.add_arg<std::string>("input.exo");
  cmdline.add_arg<std::string>("output.osh");
  cmdline.add_flag("-v", "verbose");
  if (!cmdline.parse(comm, &argc, argv) ||
      !CmdLine::check_empty(comm, argc, argv)) {
    cmdline.show_help(comm);
    return -1;
  }
  auto inpath = cmdline.get<std::string>("input.exo");
  auto outpath = cmdline.get<std::string>("output.osh");
  auto verbose = cmdline.parsed("-v");
  Omega_h::Mesh mesh(&lib);
  Omega_h::exodus::read(inpath, &mesh, verbose);
  Omega_h::binary::write(outpath, &mesh);
  return 0;
}
