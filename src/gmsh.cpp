#include "internal.hpp"

#include <fstream>

#include "algebra.hpp"
#include "classify.hpp"
#include "map.hpp"

namespace Omega_h {

namespace gmsh {

namespace {

enum {
  GMSH_VERT = 15,
  GMSH_LINE = 1,
  GMSH_TRI = 2,
  GMSH_TET = 4,
};

Int type_dim(Int type) {
  switch (type) {
    case GMSH_VERT:
      return 0;
    case GMSH_LINE:
      return 1;
    case GMSH_TRI:
      return 2;
    case GMSH_TET:
      return 3;
  }
  Omega_h_fail("omega_h can only accept linear simplices from Gmsh");
  NORETURN(-1);
}

void seek_line(std::istream& stream, std::string const& want) {
  CHECK(stream);
  std::string line;
  while (std::getline(stream, line)) {
    if (line == want) {
      break;
    }
  }
  CHECK(stream);
}

}  // end anonymous namespace

void read(std::istream& stream, Library const& lib, Mesh* mesh) {
  seek_line(stream, "$MeshFormat");
  Real format;
  Int file_type;
  Int data_size;
  stream >> format >> file_type >> data_size;
  CHECK(file_type == 0);
  CHECK(data_size == sizeof(Real));
  seek_line(stream, "$Nodes");
  LO nnodes;
  stream >> nnodes;
  CHECK(nnodes >= 0);
  std::vector<Vector<3>> node_coords;
  for (LO i = 0; i < nnodes; ++i) {
    LO number;
    stream >> number;
    // the documentation says numbers don't have to be linear,
    // but so far they have been and assuming they are saves
    // me a big lookup structure (e.g. std::map)
    CHECK(number == i + 1);
    Vector<3> coords;
    stream >> coords[0] >> coords[1] >> coords[2];
    node_coords.push_back(coords);
  }
  seek_line(stream, "$Elements");
  LO nents;
  stream >> nents;
  CHECK(nents >= 0);
  std::vector<LO> ent_class_ids[4];
  std::vector<LO> ent_nodes[4];
  for (LO i = 0; i < nents; ++i) {
    LO number;
    stream >> number;
    CHECK(number > 0);
    Int type;
    stream >> type;
    Int dim = type_dim(type);
    Int ntags;
    stream >> ntags;
    CHECK(ntags >= 2);
    Int physical;
    Int elementary;
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
  }
  CHECK(ent_nodes[2].size());
  Int max_dim = 2;
  if (ent_nodes[3].size()) {
    max_dim = 3;
  }
  HostWrite<Real> host_coords(nnodes * max_dim);
  for (LO i = 0; i < nnodes; ++i) {
    for (Int j = 0; j < max_dim; ++j) {
      host_coords[i * max_dim + j] =
          node_coords[static_cast<std::size_t>(i)][j];
    }
  }
  for (Int ent_dim = max_dim; ent_dim >= 0; --ent_dim) {
    Int neev = ent_dim + 1;
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
          mesh, lib, max_dim, eqv2v, host_coords.write());
    }
    auto eq_class_id = Read<LO>(host_class_id.write());
    LOs eq2e;
    if (ent_dim == max_dim) {
      eq2e = LOs(ndim_ents, 0, 1);
    } else if (ent_dim == VERT) {
      eq2e = eqv2v;
    } else {
      Read<I8> codes;
      auto ev2v = mesh->ask_verts_of(ent_dim);
      auto v2e = mesh->ask_up(VERT, ent_dim);
      find_matches(ent_dim, eqv2v, ev2v, v2e, &eq2e, &codes);
    }
    auto eq_class_dim = Read<I8>(ndim_ents, I8(ent_dim));
    auto class_dim =
        map_onto(eq_class_dim, eq2e, mesh->nents(ent_dim), I8(-1), 1);
    auto class_id = map_onto(eq_class_id, eq2e, mesh->nents(ent_dim), -1, 1);
    mesh->add_tag<I8>(ent_dim, "class_dim", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT, class_dim);
    mesh->add_tag<LO>(ent_dim, "class_id", 1, OMEGA_H_INHERIT, OMEGA_H_DO_OUTPUT, class_id);
  }
  project_classification(mesh);
}

void read(std::string const& filename, Library const& lib, Mesh* mesh) {
  std::ifstream file(filename.c_str());
  if (!file.is_open()) {
    Omega_h_fail("couldn't open \"%s\"\n", filename.c_str());
  }
  read(file, lib, mesh);
}
}

}  // end namespace Omega_h
