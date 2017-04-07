#ifndef OMEGA_H_ASSOC_HPP
#define OMEGA_H_ASSOC_HPP

#include <array>
#include <map>
#include <vector>
#include <string>

#include <Omega_h_defines.hpp>
#include <Omega_h_kokkos.hpp>

namespace Omega_h {

struct ClassPair {
  LO id;
  Int dim;
};

enum {
  ELEM_SET,
  SIDE_SET,
  NODE_SET,
};
enum { NSET_TYPES = 3 };

// (set_type, class_dim, set_name) -> class_ids
using GeomSets = std::array<std::array<std::map<std::string, std::vector<LO>>, DIMS>, NSET_TYPES>;
// (set_type, set_name) -> mesh_ents
using MeshSets = std::array<std::map<std::string, LOs>, NSET_TYPES>;

MeshSets invert(GeomSets const& geom_sets);

template <typename T>
OMEGA_H_DEVICE Int binary_search(Read<T> const& a, T v, LO n) {
  LO l = 0;
  LO r = n - 1;
  while (1) {
    if (l > r) return -1;
    auto m = (l + r) / 2;
    auto a_m = a[m];
    if (a_m < v) {
      l = m + 1;
      continue;
    }
    if (a_m > v) {
      r = m - 1;
      continue;
    }
    return m;
  }
}

}

#endif
