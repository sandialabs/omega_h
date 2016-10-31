#include "Omega_h.hpp"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  if (argc != 4) {
    std::cout << "usage: " << argv[0]
      << " input.osh <target nelems> output.osh\n";
    return -1;
  }
  auto path_in = argv[1];
  auto target_nelems = atoi(argv[2]);
  auto path_out = argv[3];
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(path_in, lib.world(), &mesh);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto implied_size = Omega_h::find_implied_size(&mesh);
  auto scalar = Omega_h::size_scalar_for_nelems(
      &mesh, implied_size, target_nelems);
  auto scaled_size = multiply_each_by(scalar, implied_size);
  mesh.add_tag(Omega_h::VERT, "size", 1, OMEGA_H_SIZE,
      OMEGA_H_DO_OUTPUT, scaled_size);
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;
  Omega_h::adapt(&mesh, opts);
  mesh.remove_tag(Omega_h::VERT, "size");
  Omega_h::binary::write(path_out, &mesh);
}
