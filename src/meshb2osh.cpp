#include <Omega_h.hpp>
#include <iostream>
#ifdef OMEGA_H_USE_EGADS
#include "Omega_h_egads.hpp"
#endif

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
#ifdef OMEGA_H_USE_EGADS
  if (!(argc == 3 || argc == 4)) {
    std::cout << "usage: " << argv[0]
              << " input.mesh[b] [input.egads] output.osh\n";
#else
  if (argc != 3) {
    std::cout << "usage: " << argv[0] << " input.mesh[b] output.osh\n";
#endif
    return -1;
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, argv[1]);
#ifdef OMEGA_H_USE_EGADS
  if (argc == 4) {
    auto eg = Omega_h::egads_load(argv[2]);
    Omega_h::egads_reclassify(&mesh, eg);
    Omega_h::egads_free(eg);
  }
#endif
  Omega_h::binary::write(argv[argc - 1], &mesh);
  return 0;
}
