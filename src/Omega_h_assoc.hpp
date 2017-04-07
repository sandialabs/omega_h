#ifndef OMEGA_H_ASSOC_HPP
#define OMEGA_H_ASSOC_HPP

#include <array>
#include <map>
#include <vector>
#include <string>

#include <Omega_h_defines.hpp>
#include <Omega_h_mark.hpp>

namespace Omega_h {

enum {
  ELEM_SET,
  SIDE_SET,
  NODE_SET,
};
enum { NSET_TYPES = 3 };

// (set_type, set_name) -> class_pairs
using GeomSets = std::array<std::map<std::string, std::vector<ClassPair>>, NSET_TYPES>;
// (set_type, set_name) -> mesh_ents
using MeshSets = std::array<std::map<std::string, LOs>, NSET_TYPES>;

MeshSets invert(Mesh* mesh, GeomSets const& geom_sets);

}

#endif
