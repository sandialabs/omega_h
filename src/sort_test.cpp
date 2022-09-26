#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_sort.hpp"
#include "Omega_h_file.hpp"
#include <fstream>

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
  {
    for(int i=0; i<3; i++) { 
      Read<LO> keys, gold;
      std::ifstream in("ab2b"+std::to_string(i)+".dat", std::ios::in);
      binary::read_array(in, keys, false, false);
      std::ifstream inGold("ba2ab"+std::to_string(i)+".dat", std::ios::in);
      binary::read_array(inGold, gold, false, false);
      LOs perm = sort_by_keys(keys);
      OMEGA_H_CHECK(perm == gold);
      in.close();
      inGold.close();
    }
  }
  return 0;
}
