#include "Omega_h_file.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cctype>
#include <sstream>
#include <unordered_map>

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

#ifdef OMEGA_H_USE_GMSH
#include <gmsh.h>
#endif  // OMEGA_H_USE_GMSH

namespace Omega_h {

namespace gmsh {

namespace {

enum {
  GMSH_LINE = 1,
  GMSH_TRI = 2,
  GMSH_QUAD = 3,
  GMSH_TET = 4,
  GMSH_HEX = 5,
  GMSH_VERT = 15,
};

Int type_dim(Int type) {
  switch (type) {
    case GMSH_VERT:
      return 0;
    case GMSH_LINE:
      return 1;
    case GMSH_TRI:
    case GMSH_QUAD:
      return 2;
    case GMSH_TET:
    case GMSH_HEX:
      return 3;
  }
  Omega_h_fail(
      "omega_h can only accept linear simplices and hypercubes from Gmsh");
  OMEGA_H_NORETURN(-1);
}

Omega_h_Family type_family(Int type) {
  switch (type) {
    case GMSH_VERT:
    case GMSH_LINE:
    case GMSH_TRI:
    case GMSH_TET:
      return OMEGA_H_SIMPLEX;
    case GMSH_QUAD:
    case GMSH_HEX:
      return OMEGA_H_HYPERCUBE;
  }
  OMEGA_H_NORETURN(Omega_h_Family());
}

Int gmsh_type(Omega_h_Family family, Int dim) {
  switch (family) {
    case OMEGA_H_SIMPLEX:
      switch (dim) {
        case 0:
          return GMSH_VERT;
        case 1:
          return GMSH_LINE;
        case 2:
          return GMSH_TRI;
        case 3:
          return GMSH_TET;
      }
      return -1;
    case OMEGA_H_HYPERCUBE:
      switch (dim) {
        case 0:
          return GMSH_VERT;
        case 1:
          return GMSH_LINE;
        case 2:
          return GMSH_QUAD;
        case 3:
          return GMSH_HEX;
      }
      return -1;
  }
  return -1;
}

void seek_line(std::istream& stream, std::string const& want) {
  OMEGA_H_CHECK(stream);
  std::string line;
  while (std::getline(stream, line)) {
    if (line == want) {
      break;
    }
  }
  OMEGA_H_CHECK(stream);
}

static bool seek_optional_section(
    std::istream& stream, std::string const& want) {
  OMEGA_H_CHECK(stream);
  std::string line;
  auto const pos = stream.tellg();
  bool found = false;
  while (std::getline(stream, line)) {
    if (line == want) {
      found = true;
      break;
    }
    if (!line.empty() && line.rfind("$End") != 0 && line[0] == '$') {
      // found the beginning of a new section that is not the one expected
      break;
    }
  }
  if (!found) {
    stream.clear();
    stream.seekg(pos);
  }
  OMEGA_H_CHECK(stream);
  return found;
}

static void eat_newlines(std::istream& stream) {
  while (stream.peek() == int('\n')) stream.get();
}

template <class T>
static void read(
    std::istream& stream, T& value, bool is_binary, bool needs_swapping) {
  if (is_binary) {
    binary::read_value(stream, value, needs_swapping);
  } else {
    stream >> value;
  }
}

static void read_internal_entities_section(Mesh& mesh, Real format,
    std::vector<std::string>& physical_names, std::istream& stream,
    bool is_binary, bool needs_swapping) {
  Int num_points, num_curves, num_surfaces, num_volumes;
  read(stream, num_points, is_binary, needs_swapping);
  read(stream, num_curves, is_binary, needs_swapping);
  read(stream, num_surfaces, is_binary, needs_swapping);
  read(stream, num_volumes, is_binary, needs_swapping);
  while (num_points-- > 0) {
    Int tag;
    Vector<3> point;
    Int num_physicals;
    read(stream, tag, is_binary, needs_swapping);
    read(stream, point[0], is_binary, needs_swapping);
    read(stream, point[1], is_binary, needs_swapping);
    read(stream, point[2], is_binary, needs_swapping);
    if (format == 4.0) {
      // strangely, the point is specified twice in 4.0, not 4.1
      read(stream, point[0], is_binary, needs_swapping);
      read(stream, point[1], is_binary, needs_swapping);
      read(stream, point[2], is_binary, needs_swapping);
    }
    read(stream, num_physicals, is_binary, needs_swapping);
    while (num_physicals-- > 0) {
      Int physical;
      read(stream, physical, is_binary, needs_swapping);
      OMEGA_H_CHECK(physical != 0);
      if (physical > 0) {
        const auto& physicalname = physical_names[physical - 1];
        mesh.class_sets[physicalname].emplace_back(0, tag);
      }
    }
  }
  const std::vector<std::pair<Int, Int>> params{
      {num_curves, 2}, {num_surfaces, 2}, {num_volumes, 3}};
  for (auto param : params) {
    auto num_elements = param.first;
    const auto dim = param.second;
    while (num_elements-- > 0) {
      Int tag;
      Vector<3> min_point, max_point;
      Int num_physicals;
      read(stream, tag, is_binary, needs_swapping);
      read(stream, min_point[0], is_binary, needs_swapping);
      read(stream, min_point[1], is_binary, needs_swapping);
      read(stream, min_point[2], is_binary, needs_swapping);
      read(stream, max_point[0], is_binary, needs_swapping);
      read(stream, max_point[1], is_binary, needs_swapping);
      read(stream, max_point[2], is_binary, needs_swapping);
      read(stream, num_physicals, is_binary, needs_swapping);
      while (num_physicals-- > 0) {
        Int physical;
        read(stream, physical, is_binary, needs_swapping);
        OMEGA_H_CHECK(physical != 0);
        if (physical > 0) {
          const auto& physical_name = physical_names[physical - 1];
          mesh.class_sets[physical_name].emplace_back(dim, tag);
        }
      }
      Int num_bounding_points;
      read(stream, num_bounding_points, is_binary, needs_swapping);
      while (num_bounding_points-- > 0) {
        Int points_tag;
        read(stream, points_tag, is_binary, needs_swapping);
      }
    }
  }
}

void read_internal(std::istream& stream, Mesh* mesh) {
  seek_line(stream, "$MeshFormat");
  Real format;
  Int file_type;
  Int data_size;
  stream >> format >> file_type >> data_size;
  OMEGA_H_CHECK(file_type == 0 || file_type == 1);
  bool is_binary = (file_type == 1);
  bool needs_swapping = false;
  if (is_binary) {
    eat_newlines(stream);
    int one;
    binary::read_value(stream, one, false);
    if (one != 1) {
      needs_swapping = true;
      binary::swap_bytes(one);
      OMEGA_H_CHECK(one == 1);
    }
  }
  OMEGA_H_CHECK(data_size == sizeof(Real));
  std::vector<std::string> physical_names;
  if (seek_optional_section(stream, "$PhysicalNames")) {
    Int num_physicals;
    read(stream, num_physicals, is_binary, needs_swapping);
    physical_names.reserve(static_cast<std::size_t>(num_physicals));
    eat_newlines(stream);
    for (auto i = 0; i < num_physicals; ++i) {
      Int dim, number;
      read(stream, dim, is_binary, needs_swapping);
      read(stream, number, is_binary, needs_swapping);
      OMEGA_H_CHECK(number == i + 1);
      std::string name;
      stream >> name;
      physical_names.push_back(name.substr(1, name.size() - 2));
    }
  }
  if (seek_optional_section(stream, "$Entities")) {
    read_internal_entities_section(
        *mesh, format, physical_names, stream, is_binary, needs_swapping);
    std::string line;
    std::getline(stream, line);
    // line matches "[ ]*"
    if (!line.empty()) {
      line.erase(std::remove_if(line.begin(), line.end(),
          [](unsigned char c) { return std::isspace(c); }));
    }
    OMEGA_H_CHECK(line.empty());
    std::getline(stream, line);
    OMEGA_H_CHECK(line == "$EndEntities");
  }
  seek_line(stream, "$Nodes");
  std::vector<Vector<3>> node_coords;
  std::map<int, int> node_number_map;
  int nnodes;
  if (format >= 4.0) {
    eat_newlines(stream);
    int num_entity_blocks;
    read(stream, num_entity_blocks, is_binary, needs_swapping);
    read(stream, nnodes, is_binary, needs_swapping);
    node_coords.reserve(std::size_t(nnodes));
    if (format >= 4.1) {
      int node_tag;
      read(stream, node_tag, is_binary, needs_swapping);  // min
      read(stream, node_tag, is_binary, needs_swapping);  // max
      for (int entity_block = 0; entity_block < num_entity_blocks;
           ++entity_block) {
        int class_id, class_dim;
        read(stream, class_dim, is_binary, needs_swapping);
        read(stream, class_id, is_binary, needs_swapping);
        int node_type, num_block_nodes;
        read(stream, node_type, is_binary, needs_swapping);
        read(stream, num_block_nodes, is_binary, needs_swapping);
        for (int block_node = 0; block_node < num_block_nodes; ++block_node) {
          int node_number;
          read(stream, node_number, is_binary, needs_swapping);
          const auto position = int(node_coords.size() + block_node);
          node_number_map[node_number] = position;
        }
        for (int block_node = 0; block_node < num_block_nodes; ++block_node) {
          Vector<3> coords;
          read(stream, coords[0], is_binary, needs_swapping);
          read(stream, coords[1], is_binary, needs_swapping);
          read(stream, coords[2], is_binary, needs_swapping);
          node_coords.push_back(coords);
        }
      }

    } else {
      for (int entity_block = 0; entity_block < num_entity_blocks;
           ++entity_block) {
        int class_id, class_dim;
        read(stream, class_id, is_binary, needs_swapping);
        read(stream, class_dim, is_binary, needs_swapping);
        int node_type, num_block_nodes;
        read(stream, node_type, is_binary, needs_swapping);
        read(stream, num_block_nodes, is_binary, needs_swapping);
        for (int block_node = 0; block_node < num_block_nodes; ++block_node) {
          int node_number;
          read(stream, node_number, is_binary, needs_swapping);
          node_number_map[node_number] = int(node_coords.size());
          Vector<3> coords;
          read(stream, coords[0], is_binary, needs_swapping);
          read(stream, coords[1], is_binary, needs_swapping);
          read(stream, coords[2], is_binary, needs_swapping);
          node_coords.push_back(coords);
        }
      }
    }
  } else {
    stream >> nnodes;
    OMEGA_H_CHECK(nnodes >= 0);
    node_coords.reserve(std::size_t(nnodes));
    eat_newlines(stream);
    for (LO i = 0; i < nnodes; ++i) {
      LO number;
      read(stream, number, is_binary, needs_swapping);
      // the documentation says numbers don't have to be linear,
      // but so far they have been and assuming they are saves
      // me a big lookup structure (e.g. std::map)
      OMEGA_H_CHECK(number == i + 1);
      Vector<3> coords;
      read(stream, coords[0], is_binary, needs_swapping);
      read(stream, coords[1], is_binary, needs_swapping);
      read(stream, coords[2], is_binary, needs_swapping);
      node_coords.push_back(coords);
    }
  }
  seek_line(stream, "$Elements");
  std::array<std::vector<int>, 4> ent_class_ids;
  std::array<std::vector<int>, 4> ent_nodes;
  Omega_h_Family family = OMEGA_H_SIMPLEX;
  if (format >= 4.0) {
    eat_newlines(stream);
    int num_entity_blocks, total_num_ents;
    read(stream, num_entity_blocks, is_binary, needs_swapping);
    read(stream, total_num_ents, is_binary, needs_swapping);
    if (format >= 4.1) {
      int element_tag;
      read(stream, element_tag, is_binary, needs_swapping);  // min
      read(stream, element_tag, is_binary, needs_swapping);  // max
    }
    for (int entity_block = 0; entity_block < num_entity_blocks;
         ++entity_block) {
      int class_id, class_dim;
      if (format == 4.) {
        read(stream, class_id, is_binary, needs_swapping);
        read(stream, class_dim, is_binary, needs_swapping);
      } else {
        read(stream, class_dim, is_binary, needs_swapping);
        read(stream, class_id, is_binary, needs_swapping);
      }
      int ent_type, num_block_ents;
      read(stream, ent_type, is_binary, needs_swapping);
      read(stream, num_block_ents, is_binary, needs_swapping);
      Int dim = type_dim(ent_type);
      OMEGA_H_CHECK(dim == class_dim);
      if (type_family(ent_type) == OMEGA_H_HYPERCUBE) {
        family = OMEGA_H_HYPERCUBE;
      }
      int nodes_per_ent = element_degree(family, dim, 0);
      ent_class_ids[dim].reserve(
          ent_class_ids[dim].size() + std::size_t(num_block_ents));
      ent_nodes[dim].reserve(
          ent_nodes[dim].size() + std::size_t(num_block_ents * nodes_per_ent));
      for (int block_ent = 0; block_ent < num_block_ents; ++block_ent) {
        ent_class_ids[dim].push_back(class_id);
        int ent_number;
        read(stream, ent_number, is_binary, needs_swapping);
        for (int ent_node = 0; ent_node < nodes_per_ent; ++ent_node) {
          int node_number;
          read(stream, node_number, is_binary, needs_swapping);
          auto it = node_number_map.find(node_number);
          OMEGA_H_CHECK(it != node_number_map.end());
          ent_nodes[dim].push_back(it->second);
        }
      }
    }

  } else {
    LO nents;
    stream >> nents;
    OMEGA_H_CHECK(nents >= 0);
    std::array<std::unordered_map<Int, Int>, 4> ent2physical;
    if (is_binary) {
      eat_newlines(stream);
      LO i = 0;
      while (i < nents) {
        I32 type, nfollow, ntags;
        binary::read_value(stream, type, needs_swapping);
        binary::read_value(stream, nfollow, needs_swapping);
        binary::read_value(stream, ntags, needs_swapping);
        Int dim = type_dim(type);
        if (type_family(type) == OMEGA_H_HYPERCUBE) {
          family = OMEGA_H_HYPERCUBE;
        }
        Int neev = element_degree(family, dim, 0);
        OMEGA_H_CHECK(ntags >= 2);
        for (Int j = 0; j < nfollow; ++j, ++i) {
          I32 number, physical, elementary;
          binary::read_value(stream, number, needs_swapping);
          binary::read_value(stream, physical, needs_swapping);
          binary::read_value(stream, elementary, needs_swapping);
          ent_class_ids[dim].push_back(elementary);
          if (physical != 0) {
            ent2physical[dim].emplace(elementary, physical);
          }
          for (Int k = 2; k < ntags; ++k) {
            I32 ignored;
            binary::read_value(stream, ignored, false);
          }
          for (Int k = 0; k < neev; ++k) {
            I32 node_number;
            binary::read_value(stream, node_number, needs_swapping);
            ent_nodes[dim].push_back(node_number - 1);
          }
        }
      }
    } else {
      for (LO i = 0; i < nents; ++i) {
        LO number;
        stream >> number;
        OMEGA_H_CHECK(number > 0);
        Int type;
        stream >> type;
        Int dim = type_dim(type);
        if (type_family(type) == OMEGA_H_HYPERCUBE) family = OMEGA_H_HYPERCUBE;
        Int ntags;
        stream >> ntags;
        OMEGA_H_CHECK(ntags >= 2);
        Int physical, elementary;
        stream >> physical >> elementary;
        ent_class_ids[dim].push_back(elementary);
        if (physical != 0) {
          ent2physical[dim].emplace(elementary, physical);
        }
        Int tag;
        for (Int j = 2; j < ntags; ++j) {
          stream >> tag;
        }
        Int neev = dim + 1;
        LO node_number;
        for (Int j = 0; j < neev; ++j) {
          stream >> node_number;
          ent_nodes[dim].push_back(node_number - 1);
        }
      }
    }
    for (Int dim = 0; dim < static_cast<Int>(ent2physical.size()); ++dim) {
      const auto& entities = ent2physical[dim];
      for (const auto& pair : entities) {
        const auto entity = pair.first;
        const auto physical = pair.second;
        std::string physical_name(std::to_string(physical));
        if (physical <= static_cast<Int>(physical_names.size())) {
          physical_name = physical_names[physical - 1];
        }
        mesh->class_sets[physical_name].emplace_back(dim, entity);
      }
    }
  }
  Int max_dim;
  if (ent_nodes[3].size()) {
    max_dim = 3;
  } else if (ent_nodes[2].size()) {
    max_dim = 2;
  } else if (ent_nodes[1].size()) {
    max_dim = 1;
  } else {
    Omega_h_fail("There were no Elements of dimension higher than zero!\n");
  }
  HostWrite<Real> host_coords(nnodes * max_dim);
  for (LO i = 0; i < nnodes; ++i) {
    for (Int j = 0; j < max_dim; ++j) {
      host_coords[i * max_dim + j] =
          node_coords[static_cast<std::size_t>(i)][j];
    }
  }
  for (Int ent_dim = max_dim; ent_dim >= 0; --ent_dim) {
    Int neev = element_degree(family, ent_dim, VERT);
    LO ndim_ents = static_cast<LO>(ent_nodes[ent_dim].size()) / neev;
    HostWrite<LO> host_ev2v(ndim_ents * neev);
    HostWrite<LO> host_class_id(ndim_ents);
    for (LO i = 0; i < ndim_ents; ++i) {
      for (Int j = 0; j < neev; ++j) {
        host_ev2v[i * neev + j] =
            ent_nodes[ent_dim][static_cast<std::size_t>(i * neev + j)];
      }
      host_class_id[i] = ent_class_ids[ent_dim][static_cast<std::size_t>(i)];
    }
    auto eqv2v = Read<LO>(host_ev2v.write());
    if (ent_dim == max_dim) {
      build_from_elems_and_coords(
          mesh, family, max_dim, eqv2v, host_coords.write());
    }
    classify_equal_order(mesh, ent_dim, eqv2v, host_class_id.write());
  }
  finalize_classification(mesh);
}

}  // end anonymous namespace

Mesh read(std::istream& stream, CommPtr comm) {
  auto mesh = Mesh(comm->library());
  if (comm->rank() == 0) {
    read_internal(stream, &mesh);
  }
  mesh.set_comm(comm);
  mesh.balance();
  return mesh;
}

Mesh read(filesystem::path const& filename, CommPtr comm) {
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    Omega_h_fail("couldn't open \"%s\"\n", filename.c_str());
  }
  return gmsh::read(file, comm);
}

#ifdef OMEGA_H_USE_GMSH

Mesh read_parallel(filesystem::path filename, CommPtr comm) {
  Mesh mesh(comm->library());
  Omega_h_Family family = OMEGA_H_SIMPLEX;
  {
    std::ostringstream part_suffix;
    part_suffix << '_' << comm->rank() + 1 << ".msh";
    filename += filesystem::path(part_suffix.str());
  }
  if (!exists(filename)) {
    Omega_h_fail("missing mesh part \"%s\" for rank %i\n", filename.c_str(), comm->rank());
  }
  ::gmsh::open(filename.c_str());
  std::vector<std::size_t> node_tags;
  std::vector<double> node_coords, parametric_coords;
  ::gmsh::model::mesh::getNodes(node_tags, node_coords, parametric_coords);
  const auto nnodes = static_cast<LO>(node_tags.size());
  if (nnodes == 0) {
    Omega_h_fail("Please check that filename is correct!\n");
  }

  std::map<GO, LO> node_number_map;
  Write<GO> vert_globals_w(nnodes);
  for (LO local_index = 0; local_index < nnodes; ++local_index) {
    const auto global_index =
        static_cast<GO>(node_tags[static_cast<std::size_t>(local_index)]);
    node_number_map[global_index] = local_index;
    vert_globals_w[local_index] = global_index;
  }

  std::array<std::vector<int>, 4> ent_class_ids;
  std::array<std::vector<LO>, 4> ent_nodes;
  std::array<std::vector<GO>, 4> ent_globals;
  ::gmsh::vectorpair physical_tags;
  ::gmsh::model::getPhysicalGroups(physical_tags);
  if (physical_tags.empty()) {
    physical_tags.push_back({-1, 0});
  }
  for (const auto& physical_tag : physical_tags) {
    ::gmsh::vectorpair tags;
    const auto physical = physical_tag.second;
    const auto dim = physical_tag.first;
    if (dim >= 0) {
      std::vector<int> entities_tag;
      ::gmsh::model::getEntitiesForPhysicalGroup(dim, physical, entities_tag);
      std::string name;
      ::gmsh::model::getPhysicalName(dim, physical, name);
      auto& set = mesh.class_sets[name];
      for (auto j = 0u; j < entities_tag.size(); ++j) {
        set.emplace_back(dim, entities_tag[j]);
        tags.push_back({dim, entities_tag[j]});
      }
    } else {
      ::gmsh::model::getEntities(tags);
    }
    for (auto j = 0u; j < tags.size(); ++j) {
      const auto dim = tags[j].first;
      const auto class_id = tags[j].second;
      std::vector<int> element_types;
      std::vector<std::vector<std::size_t>> element_tags, nodes_in_element_tags;
      ::gmsh::model::mesh::getElements(
          element_types, element_tags, nodes_in_element_tags, dim, class_id);
      for (std::size_t type_index = 0; type_index < element_types.size();
           ++type_index) {
        OMEGA_H_CHECK(dim == type_dim(element_types[type_index]));
        std::size_t num_block_ents = element_tags[type_index].size();
        std::size_t nodes_per_ent =
            nodes_in_element_tags[type_index].size() / num_block_ents;
        ent_class_ids[dim].reserve(ent_class_ids[dim].size() + num_block_ents);
        ent_globals[dim].reserve(ent_globals[dim].size() + num_block_ents);
        ent_nodes[dim].reserve(
            ent_nodes[dim].size() + num_block_ents * nodes_per_ent);
        for (auto block_ent = 0u; block_ent < num_block_ents; ++block_ent) {
          ent_class_ids[dim].push_back(class_id);
          ent_globals[dim].push_back(
              static_cast<GO>(element_tags[type_index][block_ent]));
          for (auto ent_node = 0u; ent_node < nodes_per_ent; ++ent_node) {
            auto it = node_number_map.find(static_cast<GO>(
                nodes_in_element_tags[type_index]
                                     [block_ent * nodes_per_ent + ent_node]));
            OMEGA_H_CHECK(it != node_number_map.end());
            ent_nodes[dim].push_back(it->second);
          }
        }
      }
    }
  }

  Int max_dim;
  if (!ent_nodes[3].empty()) {
    max_dim = 3;
  } else if (!ent_nodes[2].empty()) {
    max_dim = 2;
  } else if (!ent_nodes[1].empty()) {
    max_dim = 1;
  } else {
    Omega_h_fail("There were no Elements of dimension higher than zero!\n");
  }

  // function to decrement the values of container given in parameter
  // by the minimum value of the container across all ranks
  auto shift_container_values = [&](auto& container) {
    auto local_min_it = std::min_element(container.begin(), container.end());
    OMEGA_H_CHECK(local_min_it != container.end());
    auto global_min = comm->allreduce(*local_min_it, OMEGA_H_MIN);
    if (global_min != 0) {
      for (auto& id : container) {
        id -= global_min;
      }
    }
  };
  // shift elements global identifiers so that they start at 0
  shift_container_values(ent_globals[max_dim]);
  // shift vertices global identifiers so that they start at 0
  shift_container_values(vert_globals_w);

  HostWrite<Real> host_coords(nnodes * max_dim);
  for (LO local_index = 0; local_index < nnodes; ++local_index) {
    for (LO j = 0; j < max_dim; ++j) {
      host_coords[local_index * max_dim + j] =
          node_coords[static_cast<std::size_t>(local_index * 3 + j)];
    }
  }
  for (Int ent_dim = max_dim; ent_dim >= 0; --ent_dim) {
    const auto ndim_ents = static_cast<LO>(ent_globals[ent_dim].size());
    Write<GO> host_elem_globals(ndim_ents);
    HostWrite<LO> host_class_id(ndim_ents);
    HostWrite<LO> host_ev2v(ent_nodes[ent_dim].size());
    if (ndim_ents > 0) {
      const auto neev = static_cast<LO>(ent_nodes[ent_dim].size()) / ndim_ents;
      for (LO local_index = 0; local_index < ndim_ents; ++local_index) {
        const auto global_index =
            ent_globals[ent_dim][static_cast<std::size_t>(local_index)];
        host_elem_globals[local_index] = global_index;
        for (Int j = 0; j < neev; ++j) {
          host_ev2v[local_index * neev + j] =
              ent_nodes[ent_dim]
                       [static_cast<std::size_t>(local_index * neev + j)];
        }
        host_class_id[local_index] =
            ent_class_ids[ent_dim][static_cast<std::size_t>(local_index)];
      }
    }
    LOs eqv2v(host_ev2v.write());
    if (ent_dim == max_dim) {
      // build_from_elems2verts(
      // &mesh, mesh.library()->self(), OMEGA_H_SIMPLEX, dim, eqv2v,
      // vert_globals);
      mesh.set_comm(comm);
      mesh.set_parting(OMEGA_H_ELEM_BASED);
      mesh.set_family(family);
      mesh.set_dim(ent_dim);
      build_verts_from_globals(&mesh, vert_globals_w);
      // build_ents_from_elems2verts(&mesh, eqv2v, vert_globals, elem_globals);
      for (Int mdim = 1; mdim < ent_dim; ++mdim) {
        auto mv2v = find_unique(eqv2v, mesh.family(), ent_dim, mdim);
        add_ents2verts(&mesh, mdim, mv2v, vert_globals_w, GOs());
      }
      add_ents2verts(&mesh, ent_dim, eqv2v, vert_globals_w, host_elem_globals);
      // if (!comm->reduce_and(is_sorted(vert_globals))) {
      //  reorder_by_globals(&mesh);
      //}
      mesh.add_coords(host_coords.write());
    }
    classify_equal_order(&mesh, ent_dim, eqv2v, host_class_id.write());
  }
  mesh.set_parting(OMEGA_H_GHOSTED);
  finalize_classification(&mesh);
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  return mesh;
}

#endif  // OMEGA_H_USE_GMSH

void write(std::ostream& stream, Mesh* mesh) {
  OMEGA_H_CHECK(mesh->comm()->size() == 1);
  stream << "$MeshFormat\n";
  stream << "2.2 0 " << sizeof(Real) << '\n';
  stream << "$EndMeshFormat\n";
  stream << "$Nodes\n";
  auto nverts = mesh->nverts();
  stream << nverts << '\n';
  auto dim = mesh->dim();
  auto d_coords = mesh->coords();
  auto h_coords = HostRead<Real>(d_coords);
  stream << std::scientific << std::setprecision(17);
  for (LO i = 0; i < nverts; ++i) {
    stream << (i + 1);
    for (Int j = 0; j < dim; ++j) stream << ' ' << h_coords[i * dim + j];
    for (Int j = dim; j < 3; ++j) stream << " 0";
    stream << '\n';
  }
  stream << "$EndNodes\n";
  stream << "$Elements\n";
  LO gmsh_elem_i = 0;
  auto family = mesh->family();
  for (Int ent_dim = VERT; ent_dim <= dim; ++ent_dim) {
    auto nents = mesh->nents(ent_dim);
    auto d_class_dims = mesh->get_array<Byte>(ent_dim, "class_dim");
    auto h_class_dims = HostRead<Byte>(d_class_dims);
    for (LO i = 0; i < nents; ++i) {
      if (h_class_dims[i] == ent_dim) ++gmsh_elem_i;
    }
  }
  auto gmsh_nelems = gmsh_elem_i;
  stream << gmsh_nelems << '\n';
  gmsh_elem_i = 0;
  for (Int ent_dim = VERT; ent_dim <= dim; ++ent_dim) {
    auto nents = mesh->nents(ent_dim);
    auto d_class_ids = mesh->get_array<ClassId>(ent_dim, "class_id");
    auto d_class_dims = mesh->get_array<Byte>(ent_dim, "class_dim");
    auto d_ents2verts = mesh->ask_verts_of(ent_dim);
    auto h_class_ids = HostRead<ClassId>(d_class_ids);
    auto h_class_dims = HostRead<Byte>(d_class_dims);
    auto h_ents2verts = HostRead<LO>(d_ents2verts);
    auto deg = element_degree(family, ent_dim, VERT);
    for (LO i = 0; i < nents; ++i) {
      if (h_class_dims[i] != ent_dim) continue;
      stream << ++gmsh_elem_i;
      stream << ' ' << gmsh_type(family, ent_dim);
      stream << " 2";
      auto class_id = h_class_ids[i];
      stream << ' ' << class_id << ' ' << class_id;
      for (Int j = 0; j < deg; ++j) {
        stream << ' ' << (h_ents2verts[i * deg + j] + 1);
      }
      stream << '\n';
    }
  }
  stream << "$EndElements\n";
}

void write(filesystem::path const& filepath, Mesh* mesh) {
  std::ofstream stream(filepath.c_str());
  gmsh::write(stream, mesh);
}

#ifdef OMEGA_H_USE_GMSH

void write_parallel(filesystem::path const& filename, Mesh& mesh) {
  const auto dim = mesh.dim();
  const auto nnodes = mesh.nverts();
  const auto globals_v = mesh.globals(VERT);
  const auto coords = mesh.coords();
  const auto gmsh_dims = 3;

  std::vector<std::size_t> node_tags(static_cast<size_t>(nnodes));
  std::vector<double> node_coords(gmsh_dims * static_cast<size_t>(nnodes), 0);
  for (LO local_index = 0; local_index < nnodes; ++local_index) {
    node_tags[static_cast<std::size_t>(local_index)] =
        static_cast<size_t>(globals_v[local_index]) + 1;
    // cleaner to use gather_vector on coords
    for (LO j = 0; j < dim; ++j) {
      node_coords[static_cast<std::size_t>(local_index * gmsh_dims + j)] =
          coords[local_index * dim + j];
    }
  }

  auto ents2verts = mesh.ask_elem_verts();
  const LO ndim_ents = mesh.nelems();
  std::vector<std::size_t> element_tags(static_cast<size_t>(ndim_ents));
  std::vector<std::size_t> ent_nodes(
      static_cast<size_t>((dim + 1) * ndim_ents));
  const auto globals_e = mesh.globals(dim);
  for (LO local_index = 0; local_index < ndim_ents; ++local_index) {
    element_tags[static_cast<std::size_t>(local_index)] =
        static_cast<size_t>(globals_e[local_index]) + 1;
    for (Int j = 0; j < dim + 1; ++j) {
      ent_nodes[static_cast<std::size_t>(local_index * (dim + 1) + j)] =
          static_cast<size_t>(
              globals_v[ents2verts[local_index * (dim + 1) + j]]) +
          1;
    }
  }
  ::gmsh::model::add("parallel_export");

  // add discrete geometry with tag 1
  ::gmsh::model::addDiscreteEntity(dim, 1);
  ::gmsh::model::mesh::addNodes(dim, 1, node_tags, node_coords);
  ::gmsh::model::mesh::addElementsByType(
      1, dim == 2 ? 2 : 4, element_tags, ent_nodes);
  std::ostringstream partFilename;
  partFilename << filename << '_' << mesh.comm()->rank() + 1 << ".msh";
  ::gmsh::write(partFilename.str());
}

#endif  // OMEGA_H_USE_GMSH

}  // namespace gmsh

}  // end namespace Omega_h
