#include "Omega_h_library.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_sort.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_atomics.hpp"
#include "Omega_h_for.hpp"
#include <fstream>

int main(int argc, char** argv) {
  using namespace Omega_h;
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  {
    LOs a({0, 2, 0, 1});
    LOs perm = sort_by_keys(a,1);
    LOs gold({0, 2, 3, 1});
    auto perm_hr = HostRead<LO>(perm);
    auto gold_hr = HostRead<LO>(gold);
    for(int j=0; j<perm_hr.size(); j++) {
      fprintf(stderr, "%d %d %d\n", j, perm_hr[j], gold_hr[j]);
    }
    OMEGA_H_CHECK(perm == gold);
  }
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
      fprintf(stderr, "large test %d\n", i);
      Read<LO> keys, gold;
      std::ifstream in("ab2b"+std::to_string(i)+".dat", std::ios::in);
      assert(in.is_open());
      binary::read_array(in, keys, false, false);
      std::ifstream inGold("ba2ab"+std::to_string(i)+".dat", std::ios::in);
      assert(in.is_open());
      binary::read_array(inGold, gold, false, false);
      in.close();
      inGold.close();
      LOs perm = sort_by_keys(keys);
      auto perm_hr = HostRead<LO>(perm);
      auto gold_hr = HostRead<LO>(gold);
      bool isSame = true;
      assert(perm_hr.size() == gold_hr.size());
      for(int j=0; j<perm_hr.size(); j++) {
        if(perm_hr[j] != gold_hr[j]) {
          isSame = false;
          fprintf(stderr, "%d %d %d\n", j, perm_hr[j], gold_hr[j]);
        }
      }
      fprintf(stderr, "host matches %s\n", (isSame) ? "yes" : "no");
      Write<LO> cnt({0});
      auto countNEQ = OMEGA_H_LAMBDA(int i) {
        if(perm[i] != gold[i]) {
          atomic_increment(&cnt[0]);
        }
      };
      parallel_for(perm.size(), countNEQ);
      auto cnt_hr = HostRead<LO>(cnt);
      fprintf(stderr, "device matches %s\n", (cnt_hr[0] == 0) ? "yes" : "no");
      auto permMatch = (perm == gold);
      fprintf(stderr, "perm matches (==) %s\n", (permMatch) ? "yes" : "no");
      OMEGA_H_CHECK(permMatch);
    }
  }
  return 0;
}
