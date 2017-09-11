#include "Omega_h.hpp"
#include "Omega_h_array_ops.hpp"

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
  auto target_nelems = atof(argv[2]);
  auto path_out = argv[3];
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(path_in, lib.world(), &mesh);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto metrics = Omega_h::get_implied_isos(&mesh);
  auto scalar =
      Omega_h::get_metric_scalar_for_nelems(&mesh, metrics, target_nelems);
  metrics = multiply_each_by(metrics, scalar);
  mesh.add_tag(Omega_h::VERT, "metric", 1, metrics);
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;
  Omega_h::adapt(&mesh, opts);
  mesh.remove_tag(Omega_h::VERT, "metric");
  Omega_h::binary::write(path_out, &mesh);
}
