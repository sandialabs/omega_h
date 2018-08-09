#include <cstdlib>
#include <iostream>

#include <Omega_h_build.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_for.hpp>
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
  auto& factor_flag = cmdline.add_flag("--factor", "weight for saltzman shift");
  factor_flag.add_arg<double>("factor");
  if (!cmdline.parse_final(world, &argc, argv)) return -1;
  auto x = cmdline.get<double>("length");
  auto y = cmdline.get<double>("width");
  auto z = cmdline.get<double>("height");
  auto nx = cmdline.get<int>("nx");
  auto ny = cmdline.get<int>("ny");
  auto nz = cmdline.get<int>("nz");
  OMEGA_H_CHECK(nz > 0);
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
  auto coords = mesh.coords();
  auto dx = x / nx;
  auto pi_over_x = Omega_h::PI / x;
  Omega_h::Real factor = 1.;
  if (cmdline.parsed("--factor"))
    factor = cmdline.get<double>("--factor", "factor");
  auto new_coords = Omega_h::Write<Omega_h::Real>(mesh.nverts() * 3);
  auto f = OMEGA_H_LAMBDA(Omega_h::LO v) {
    auto v_pos = Omega_h::get_vector<3>(coords, v);
    auto v_x = v_pos[0];
    auto v_y = v_pos[1];
    auto v_z = v_pos[2];
    auto v_y_norm = (v_y / y) - (1. / 2.);
    auto v_z_norm = (v_z / z) - (1. / 2.);
    auto v_yz_norm = v_y_norm * v_z_norm;
    auto x_shift = factor * v_yz_norm * dx * std::sin(pi_over_x * v_x);
    v_pos[0] += x_shift;
    Omega_h::set_vector(new_coords, v, v_pos);
  };
  Omega_h::parallel_for(mesh.nverts(), f);
  mesh.set_coords(Omega_h::Reals(new_coords));
  Omega_h::binary::write(outdir, &mesh);
  return 0;
}
