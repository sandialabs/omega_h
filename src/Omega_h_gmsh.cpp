#include "Omega_h_file.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>

#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_vector.hpp"

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

static void register_physical_entity(Mesh& mesh,
    const std::vector<std::string>& names, Int dim, Int number, Int physical) {
  if (physical != 0) {
    auto& set =
        mesh.class_sets[names.at(static_cast<std::size_t>(physical - 1))];
    set.emplace_back(dim, number - 1);
  }
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

static void read_internal_entities_section(Real format,
    std::array<std::map<Int, std::vector<Int>>, 4>& entity2physicals,
    std::istream& stream, bool is_binary, bool needs_swapping) {
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
      entity2physicals[0][tag].push_back(physical);
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
        entity2physicals[dim][tag].push_back(physical);
      }
      Int num_bounding_points;
      read(stream, num_bounding_points, is_binary, needs_swapping);
      while (num_bounding_points-- > 0) {
        Int tag;
        read(stream, tag, is_binary, needs_swapping);
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
  std::array<std::map<Int, std::vector<Int>>, 4> entity2physicals;
  if (seek_optional_section(stream, "$Entities")) {
    read_internal_entities_section(
        format, entity2physicals, stream, is_binary, needs_swapping);
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
  std::vector<int> ent_class_ids[4];
  std::vector<int> ent_nodes[4];
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
        ent_nodes[dim].reserve(ent_nodes[dim].size() +
                               std::size_t(num_block_ents * nodes_per_ent));
        for (int block_ent = 0; block_ent < num_block_ents; ++block_ent) {
          ent_class_ids[dim].push_back(class_id);
          int ent_number;
          read(stream, ent_number, is_binary, needs_swapping);
          if (class_id != 0) {
            for (const auto physical : entity2physicals[class_dim][class_id]) {
              register_physical_entity(
                  *mesh, physical_names, class_dim, ent_number, physical);
            }
          }
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
          for (Int k = 2; k < ntags; ++k) {
            I32 ignored;
            binary::read_value(stream, ignored, false);
          }
          for (Int k = 0; k < neev; ++k) {
            I32 node_number;
            binary::read_value(stream, node_number, needs_swapping);
            ent_nodes[dim].push_back(node_number - 1);
          }
          register_physical_entity(
              *mesh, physical_names, dim, number, physical);
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
        register_physical_entity(*mesh, physical_names, dim, number, physical);
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

}  // namespace gmsh

}  // end namespace Omega_h
