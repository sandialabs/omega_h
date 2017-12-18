#include <iostream>

#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  if (argc != 5) {
    std::cout << "usage: " << argv[0]
              << " input.osh input.solb <sol-name> output.osh\n";
    return -1;
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], lib.world(), &mesh);
  Omega_h::meshb::read_sol(&mesh, argv[2], argv[3]);
  Omega_h::binary::write(argv[4], &mesh);
  return 0;
}
