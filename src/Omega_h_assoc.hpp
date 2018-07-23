#ifndef OMEGA_H_ASSOC_HPP
#define OMEGA_H_ASSOC_HPP

#include <array>
#include <string>
#include <vector>

#include <Omega_h_defines.hpp>
#include <Omega_h_mark.hpp>

namespace Omega_h {

enum SetType {
  ELEM_SET,
  SIDE_SET,
  NODE_SET,
};
enum { NSET_TYPES = 3 };

extern char const* const assoc_names[NSET_TYPES];

// (set_type, set_name) -> class_pairs
using Assoc = std::array<ClassSets, NSET_TYPES>;
using MeshDimSets = std::map<std::string, LOs>;
// (set_type, set_name) -> mesh_ents
using MeshSets = std::array<MeshDimSets, NSET_TYPES>;

Int get_assoc_dim(size_t set_type, Int mesh_dim);

MeshSets invert(Mesh* mesh, Assoc const& geom_sets);

void update_from_file(Assoc* p_assoc, std::string const& filepath);

Assoc read_assoc_file(std::string const& filepath);

Assoc get_box_assoc(int dim);

}  // namespace Omega_h

#endif
