#include <cstdlib>
#include <iostream>

#include <Omega_h_library.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  if (argc != 8) {
    if (world->rank() == 0) {
      std::cout << "usage: " << argv[0]
                << " length width height nx ny nz output.osh\n";
      std::cout << " where ny is the number of elements along the Y axis.\n";
      std::cout << " set nz=0 to generate a 2D mesh.\n";
    }
    return -1;
  }
  auto x = atof(argv[1]);
  auto y = atof(argv[2]);
  auto z = atof(argv[3]);
  auto nx = atoi(argv[4]);
  auto ny = atoi(argv[5]);
  auto nz = atoi(argv[6]);
  auto outdir = argv[7];
  auto mesh = Omega_h::build_box(world, x, y, z, nx, ny, nz);
  Omega_h::binary::write(outdir, &mesh);
  return 0;
}
