#include <Omega_h_assoc.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>

namespace Omega_h {

static Int set_ent_dim(size_t set_type, Int mesh_dim) {
  switch (set_type) {
    case ELEM_SET: return mesh_dim;
    case SIDE_SET: return mesh_dim - 1;
    case NODE_SET: return VERT;
  }
  return -1;
}

MeshSets invert(Mesh* mesh, GeomSets const& geom_sets) {
  MeshSets mesh_sets;
  for (size_t set_type = 0; set_type < NSET_TYPES; ++set_type) {
    auto ent_dim = set_ent_dim(set_type, mesh->dim());
    for (auto& name_pairs : geom_sets[set_type]) {
      auto& name = name_pairs.first;
      auto& pairs = name_pairs.second;
      auto marks = mark_class_closures(mesh, ent_dim, pairs);
      auto set_ents = collect_marked(marks);
      mesh_sets[set_type][name] = set_ents;
    }
  }
  return mesh_sets;
}

}
