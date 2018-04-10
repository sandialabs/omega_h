#include <cstdlib>
#include <iostream>

#include <Omega_h_build.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<double>("length");
  cmdline.add_arg<double>("width");
  cmdline.add_arg<double>("height");
  cmdline.add_arg<int>("nx");
  cmdline.add_arg<int>("ny");
  cmdline.add_arg<int>("nz");
  cmdline.add_arg<std::string>("output.osh");
  auto& family_flag = cmdline.add_flag("--family", "simplex or hypercube");
  cmdline.add_flag("--symmetric", "split hypercubes symmetrically");
  family_flag.add_arg<std::string>("type");
  if (!cmdline.parse_final(world, &argc, argv)) return -1;
  auto x = cmdline.get<double>("length");
  auto y = cmdline.get<double>("width");
  auto z = cmdline.get<double>("height");
  auto nx = cmdline.get<int>("nx");
  auto ny = cmdline.get<int>("ny");
  auto nz = cmdline.get<int>("nz");
  auto outdir = cmdline.get<std::string>("output.osh");
  auto family = OMEGA_H_SIMPLEX;
  if (cmdline.parsed("--family")) {
    auto family_type = cmdline.get<std::string>("--family", "type");
    if (family_type == "simplex")
      family = OMEGA_H_SIMPLEX;
    else if (family_type == "hypercube")
      family = OMEGA_H_HYPERCUBE;
    else {
      std::cout << "unknown family: " << family_type << std::endl;
      cmdline.show_help(world, argv);
      return -1;
    }
  }
  auto symmetric = cmdline.parsed("--symmetric");
  auto mesh = Omega_h::build_box(world, family, x, y, z, nx, ny, nz, symmetric);
  Omega_h::binary::write(outdir, &mesh);
  return 0;
}
