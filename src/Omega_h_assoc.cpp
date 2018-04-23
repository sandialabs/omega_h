#include <Omega_h_assoc.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mesh.hpp>

#include <fstream>
#include <sstream>

namespace Omega_h {

static char const* const assoc_file_names[NSET_TYPES] = {
    "element set", "side set", "node set"};

Int get_assoc_dim(size_t set_type, Int mesh_dim) {
  switch (set_type) {
    case ELEM_SET:
      return mesh_dim;
    case SIDE_SET:
      return mesh_dim - 1;
    case NODE_SET:
      return VERT;
  }
  return -1;
}

MeshSets invert(Mesh* mesh, Assoc const& geom_sets) {
  MeshSets mesh_sets;
  for (size_t set_type = 0; set_type < NSET_TYPES; ++set_type) {
    auto ent_dim = get_assoc_dim(set_type, mesh->dim());
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

void update_from_file(Assoc* p_assoc, std::string const& filepath) {
  Assoc& assoc = *p_assoc;
  std::ifstream f(filepath.c_str());
  if (!f.is_open()) {
    Omega_h_fail("Could not open associations file \"%s\"\n", filepath.c_str());
  }
  std::string sline;
  LO lc = 0;
  while (std::getline(f, sline)) {
    if (sline.empty()) break;
    ++lc;
    std::string rest;
    size_t set_type;
    for (set_type = 0; set_type < NSET_TYPES; ++set_type) {
      std::string set_name = assoc_file_names[set_type];
      if (sline.compare(0, set_name.length(), set_name) == 0) {
        rest = sline.substr(set_name.length());
        break;
      }
    }
    if (set_type >= NSET_TYPES) {
      Omega_h_fail("Unknown set type \"%s\" at %s +%d\n", sline.c_str(),
          filepath.c_str(), lc);
    }
    std::stringstream rest_stream(rest);
    std::string set_name;
    rest_stream >> set_name;
    LO set_size;
    rest_stream >> set_size;
    if (!rest_stream) {
      Omega_h_fail(
          "Couldn't parse set name and size at %s +%d\n", filepath.c_str(), lc);
    }
    for (LO i = 0; i < set_size; ++i) {
      std::string eline;
      std::getline(f, eline);
      if (!f || eline.empty()) {
        Omega_h_fail(
            "Expected more pairs after %s +%d\n", filepath.c_str(), lc);
      }
      ++lc;
      std::stringstream pair_stream(eline);
      Int class_dim;
      LO class_id;
      pair_stream >> class_dim >> class_id;
      if (!pair_stream) {
        Omega_h_fail("Couldn't parse pair \"%s\" at %s +%d\n", eline.c_str(),
            filepath.c_str(), lc);
      }
      assoc[set_type][set_name].push_back({class_dim, class_id});
    }
  }
}

Assoc read_assoc_file(std::string const& filepath) {
  Assoc assoc;
  update_from_file(&assoc, filepath);
  return assoc;
}

Assoc get_box_assoc(int dim) {
  Assoc assoc;
  if (dim == 1) {
    assoc[NODE_SET]["x-"] = {{0, 0}};
    assoc[SIDE_SET]["x-"] = {{0, 0}};
    assoc[NODE_SET]["body"] = {{1, 1}};
    assoc[ELEM_SET]["body"] = {{1, 1}};
    assoc[NODE_SET]["x+"] = {{0, 2}};
    assoc[SIDE_SET]["x+"] = {{0, 2}};
  } else if (dim == 2) {
    assoc[NODE_SET]["y-"] = {{1, 1}};
    assoc[SIDE_SET]["y-"] = {{1, 1}};
    assoc[NODE_SET]["x-"] = {{1, 3}};
    assoc[SIDE_SET]["x-"] = {{1, 3}};
    assoc[NODE_SET]["body"] = {{2, 4}};
    assoc[ELEM_SET]["body"] = {{2, 4}};
    assoc[NODE_SET]["x+"] = {{1, 5}};
    assoc[SIDE_SET]["x+"] = {{1, 5}};
    assoc[NODE_SET]["y+"] = {{1, 7}};
    assoc[SIDE_SET]["y+"] = {{1, 7}};
  } else if (dim == 3) {
    assoc[NODE_SET]["z-"] = {{2, 4}};
    assoc[SIDE_SET]["z-"] = {{2, 4}};
    assoc[NODE_SET]["y-"] = {{2, 10}};
    assoc[SIDE_SET]["y-"] = {{2, 10}};
    assoc[NODE_SET]["x-"] = {{2, 12}};
    assoc[SIDE_SET]["x-"] = {{2, 12}};
    assoc[NODE_SET]["body"] = {{3, 13}};
    assoc[ELEM_SET]["body"] = {{3, 13}};
    assoc[NODE_SET]["x+"] = {{2, 14}};
    assoc[SIDE_SET]["x+"] = {{2, 14}};
    assoc[NODE_SET]["y+"] = {{2, 16}};
    assoc[SIDE_SET]["y+"] = {{2, 16}};
    assoc[NODE_SET]["z+"] = {{2, 22}};
    assoc[SIDE_SET]["z+"] = {{2, 22}};
  }
  return assoc;
}

}  // namespace Omega_h
