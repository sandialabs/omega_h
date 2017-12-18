#include <iostream>

#include <Omega_h_library.hpp>
#include <Omega_h_timer.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  if (argc != 4) {
    if (!world->rank()) {
      std::cout << "usage: " << argv[0] << " in.osh <nparts> out.osh\n";
    }
    return -1;
  }
  auto nparts_total = world->size();
  auto path_in = argv[1];
  auto nparts_out = atoi(argv[2]);
  auto path_out = argv[3];
  auto t0 = Omega_h::now();
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
  auto nparts_in = Omega_h::binary::read_nparts(path_in, world);
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
  auto mesh = Omega_h::Mesh(&lib);
  auto version = Omega_h::binary::read_version(path_in, world);
  if (is_in) {
    Omega_h::binary::read_in_comm(path_in, comm_in, &mesh, version);
    if (nparts_out < nparts_in) {
      Omega_h_fail(
          "partitioning to a smaller part count not yet implemented\n");
    }
  }
  if (is_in || is_out) mesh.set_comm(comm_out);
  if (is_out) {
    if (nparts_out != nparts_in) mesh.balance();
    Omega_h::binary::write(path_out, &mesh);
  }
  world->barrier();
  auto t1 = Omega_h::now();
  auto imb = mesh.imbalance();
  if (!world->rank()) {
    std::cout << "repartitioning took " << (t1 - t0) << " seconds\n";
    std::cout << "imbalance is " << imb << "\n";
  }
}
