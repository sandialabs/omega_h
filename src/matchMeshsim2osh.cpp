#include <iostream>

#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();

  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh-in");
  cmdline.add_arg<std::string>("model-in");
  cmdline.add_arg<std::string>("mesh-out");
  cmdline.add_arg<int>("nparts_out");

  auto nparts_total = world->size();
  auto nparts_in = 1;//input distributed mesh not implemented

  if (!cmdline.parse_final(world, &argc, argv)) return -1;
  auto mesh_in = cmdline.get<std::string>("matched-simmodsuite-mesh-in");
  auto model_in = cmdline.get<std::string>("simmodsuite-geomsim-model-in");
  auto mesh_out = cmdline.get<std::string>("mesh-out");
  auto nparts_out = cmdline.get<int>("nparts_out");
  auto mesh = Omega_h::Mesh(&lib);

  if (nparts_out < 1) {
    if (!world->rank()) {
      std::cout << "error: invalid output part count " << nparts_out << '\n';
    }
    return -1;
  }
  if (nparts_out > nparts_total) {
    if (!world->rank()) {
      std::cout << "error: output part count " << nparts_out
                << " greater than MPI job size " << nparts_total << '\n';
    }
    return -1;
  }
  if (nparts_in > nparts_total) {
    if (!world->rank()) {
      std::cout << "error: input part count " << nparts_in
                << " greater than MPI job size " << nparts_total << '\n';
    }
    return -1;
  }

  auto is_in = (world->rank() < nparts_in);
  auto comm_in = world->split(int(is_in), 0);
  auto is_out = (world->rank() < nparts_out);
  auto comm_out = world->split(int(is_out), 0);

  Omega_h::meshsim::matchRead(mesh_in, model_in, comm_in, &mesh, is_in);
  if (is_in || is_out) {
    mesh.set_comm(comm_out);
  }
  if (is_out) {
    if (nparts_out != nparts_in) mesh.balance();

    Omega_h::binary::write(mesh_out, &mesh);
  }

  return 0;
}
