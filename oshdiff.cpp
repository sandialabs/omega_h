#include "omega_h.hpp"

#include <iostream>

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  Real tol = 1e-6;
  Real floor = 0.0;
  bool get_tol = false;
  bool get_floor = false;
  bool get_help = false;
  char const* filea = nullptr;
  char const* fileb = nullptr;
  bool allow_superset = false;
  for (int i = 1; i < argc; ++i) {
    if (get_tol) {
      tol = atof(argv[i]);
      get_tol = false;
      continue;
    }
    if (get_floor) {
      floor = atof(argv[i]);
      get_floor = false;
      continue;
    }
    if (!strcmp(argv[i], "-tolerance")) {
      get_tol = true;
      continue;
    }
    if (!strcmp(argv[i], "-Floor")) {
      get_floor = true;
      continue;
    }
    if (!strcmp(argv[i], "-help")) {
      get_help = true;
      continue;
    }
    if (!strcmp(argv[i], "-superset")) {
      allow_superset = true;
      continue;
    }
    if (!filea) {
      filea = argv[i];
      continue;
    }
    if (!fileb) {
      fileb = argv[i];
      continue;
    }
  }
  if (get_tol) {
    printf("-tolerance needs an argument\n");
    return -1;
  }
  if (get_floor) {
    printf("-Floor needs an argument\n");
    return -1;
  }
  if (get_help || !filea || !fileb) {
    std::cout << "\n";
    std::cout << "usage: oshdiff [options] dir1.osh dir2.osh\n";
    std::cout << "   or: oshdiff [-help]             (usage)\n";
    std::cout << "\n";
    std::cout << "    -help (Print this summary and exit.)\n";
    std::cout << "    -tolerance <$val> (Overrides the default tolerance of 1.0E-6.)\n";
    std::cout << "    -Floor <$val> (Overrides the default floor tolerance of 0.0.)\n";
    return -1;
  }
  osh::Mesh a;
  osh::binary::read(filea, lib.world(), &a);
  osh::Mesh b;
  osh::binary::read(fileb, lib.world(), &b);
  auto res = compare_meshes(&a, &b, tol, floor, true);
  if (res == SAME_MESH) return 0;
  if (allow_superset && res == SUPERSET_MESH) return 0;
  return 2;
}
