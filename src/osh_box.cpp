#include <cstdlib>
#include <iostream>

#include "Omega_h.hpp"
#include "box.hpp"

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
  Omega_h::Mesh mesh(&lib);
  if (world->rank() == 0) {
    Omega_h::build_box(&mesh, x, y, z, nx, ny, nz);
    Omega_h::classify_by_angles(&mesh, Omega_h::PI / 4);
    Omega_h::set_box_class_ids(&mesh, x, y, z, nx, ny, nz);
  }
  mesh.set_comm(world);
  mesh.balance();
  Omega_h::binary::write(outdir, &mesh);
  return 0;
}
