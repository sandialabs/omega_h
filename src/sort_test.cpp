#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_sort.hpp"

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  {
    LOs a({0, 2, 0, 1});
    LOs perm = sort_by_keys(a, 2);
    OMEGA_H_CHECK(perm == LOs({1, 0}));
  }
  {
    LOs a({0, 2, 1, 1});
    LOs perm = sort_by_keys(a, 2);
    OMEGA_H_CHECK(perm == LOs({0, 1}));
  }
  {
    LOs a({1, 2, 3, 1, 2, 2, 3, 0, 0});
    LOs perm = sort_by_keys(a, 3);
    OMEGA_H_CHECK(perm == LOs({1, 0, 2}));
  }
  return 0;
}
