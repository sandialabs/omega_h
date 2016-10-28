#include <Omega_h.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  if (argc != 3) {
    std::cout << "usage: " << argv[0] << " input.meshb output.osh\n";
    return -1;
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, argv[1]);
  Omega_h::binary::write(argv[2], &mesh);
  return 0;
}
