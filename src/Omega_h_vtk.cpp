#include "Omega_h_vtk.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>

#ifdef OMEGA_H_USE_ZLIB
#include <zlib.h>
#endif

#include "Omega_h_base64.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_simplex.hpp"
#include "Omega_h_tag.hpp"
#include "Omega_h_xml.hpp"

namespace Omega_h {

namespace vtk {

TagSet get_all_vtk_tags(Mesh* mesh) {
  auto out = get_all_mesh_tags(mesh);
  for (Int i = 0; i < mesh->dim(); ++i) {
    out[size_t(i)].insert("local");
    if (mesh->comm()->size() > 1) {
      out[size_t(i)].insert("owner");
    }
  }
  return out;
}

namespace {

/* start of C++ ritual dance to print a string based on
   type properties */

template <bool is_signed, std::size_t size>
struct IntTraits;

template <>
struct IntTraits<true, 1> {
  inline static char const* name() { return "Int8"; }
};

template <>
struct IntTraits<true, 4> {
  inline static char const* name() { return "Int32"; }
};

template <>
struct IntTraits<true, 8> {
  inline static char const* name() { return "Int64"; }
};

template <>
struct IntTraits<false, 8> {
  inline static char const* name() { return "UInt64"; }
};

template <std::size_t size>
struct FloatTraits;

template <>
struct FloatTraits<8> {
  inline static char const* name() { return "Float64"; }
};

template <typename T, typename Enable = void>
struct Traits;

template <typename T>
struct Traits<T, typename std::enable_if<std::is_integral<T>::value>::type> {
  inline static char const* name() {
    return IntTraits<std::is_signed<T>::value, sizeof(T)>::name();
  }
};

template <typename T>
struct Traits<T,
    typename std::enable_if<std::is_floating_point<T>::value>::type> {
  inline static char const* name() { return FloatTraits<sizeof(T)>::name(); }
};

/* end of C++ ritual dance to get a string based on type properties */

template <typename T>
void describe_array(std::ostream& stream, std::string const& name, Int ncomps) {
  stream << "type=\"" << Traits<T>::name() << "\"";
  stream << " Name=\"" << name << "\"";
  stream << " NumberOfComponents=\"" << ncomps << "\"";
  stream << " format=\"binary\"";
}

bool read_array_start_tag(std::istream& stream, Omega_h_Type* type_out,
    std::string* name_out, Int* ncomps_out) {
  auto st = xml::read_tag(stream);
  if (st.elem_name != "DataArray" || st.type != xml::Tag::START) {
    OMEGA_H_CHECK(st.type == xml::Tag::END);
    return false;
  }
  auto type_name = st.attribs["type"];
  if (type_name == "Int8")
    *type_out = OMEGA_H_I8;
  else if (type_name == "Int32")
    *type_out = OMEGA_H_I32;
  else if (type_name == "Int64")
    *type_out = OMEGA_H_I64;
  else if (type_name == "Float64")
    *type_out = OMEGA_H_F64;
  *name_out = st.attribs["Name"];
  *ncomps_out = std::stoi(st.attribs["NumberOfComponents"]);
  OMEGA_H_CHECK(st.attribs["format"] == "binary");
  return true;
}

template <typename T>
void write_array(
    std::ostream& stream, std::string const& name, Int ncomps, Read<T> array) {
  if (!(array.exists())) {
    Omega_h_fail("vtk::write_array: \"%s\" doesn't exist\n", name.c_str());
  }
  stream << "<DataArray ";
  describe_array<T>(stream, name, ncomps);
  stream << ">\n";
  HostRead<T> uncompressed(array);
  std::uint64_t uncompressed_bytes =
      sizeof(T) * static_cast<uint64_t>(array.size());
#ifdef OMEGA_H_USE_ZLIB
  uLong source_bytes = uncompressed_bytes;
  uLong dest_bytes = ::compressBound(source_bytes);
  auto compressed = new ::Bytef[dest_bytes];
  int ret = ::compress2(compressed, &dest_bytes,
      reinterpret_cast<const ::Bytef*>(nonnull(uncompressed.data())),
      source_bytes, Z_BEST_SPEED);
  OMEGA_H_CHECK(ret == Z_OK);
  std::string encoded = base64::encode(compressed, dest_bytes);
  delete[] compressed;
  std::uint64_t header[4] = {
      1, uncompressed_bytes, uncompressed_bytes, dest_bytes};
  std::string enc_header = base64::encode(header, sizeof(header));
#else
  std::string enc_header =
      base64::encode(&uncompressed_bytes, sizeof(std::uint64_t));
  std::string encoded =
      base64::encode(nonnull(uncompressed.data()), uncompressed_bytes);
#endif
  stream << enc_header << encoded << '\n';
  stream << "</DataArray>\n";
}

template <typename T>
Read<T> read_array(
    std::istream& stream, LO size, bool is_little_endian, bool is_compressed) {
  auto enc_both = base64::read_encoded(stream);
  std::uint64_t uncompressed_bytes, compressed_bytes;
  std::string encoded;
#ifdef OMEGA_H_USE_ZLIB
  if (is_compressed) {
    std::uint64_t header[4];
    auto nheader_chars = base64::encoded_size(sizeof(header));
    auto enc_header = enc_both.substr(0, nheader_chars);
    base64::decode(enc_header, header, sizeof(header));
    for (std::uint64_t i = 0; i < 4; ++i) {
      binary::swap_if_needed(header[i], is_little_endian);
    }
    encoded = enc_both.substr(nheader_chars);
    uncompressed_bytes = header[2];
    compressed_bytes = header[3];
  } else
#else
  OMEGA_H_CHECK(is_compressed == false);
#endif
  {
    auto nheader_chars = base64::encoded_size(sizeof(std::uint64_t));
    auto enc_header = enc_both.substr(0, nheader_chars);
    base64::decode(enc_header, &uncompressed_bytes, sizeof(uncompressed_bytes));
    binary::swap_if_needed(uncompressed_bytes, is_little_endian);
    compressed_bytes = uncompressed_bytes;
    encoded = enc_both.substr(nheader_chars);
  }
  OMEGA_H_CHECK(uncompressed_bytes == std::uint64_t(size) * sizeof(T));
  HostWrite<T> uncompressed(size);
#ifdef OMEGA_H_USE_ZLIB
  if (is_compressed) {
    auto compressed = new ::Bytef[compressed_bytes];
    base64::decode(encoded, compressed, compressed_bytes);
    uLong dest_bytes = static_cast<uLong>(uncompressed_bytes);
    uLong source_bytes = static_cast<uLong>(compressed_bytes);
    ::Bytef* uncompressed_ptr =
        reinterpret_cast< ::Bytef*>(nonnull(uncompressed.data()));
    int ret =
        ::uncompress(uncompressed_ptr, &dest_bytes, compressed, source_bytes);
    if (ret != Z_OK) {
      Omega_h_fail(
          "code %d: couln't decompress \"%s\"\n", ret, encoded.c_str());
    }
    OMEGA_H_CHECK(dest_bytes == static_cast<uLong>(uncompressed_bytes));
    delete[] compressed;
  } else
#endif
  {
    base64::decode(encoded, nonnull(uncompressed.data()), uncompressed_bytes);
  }
  return binary::swap_if_needed(
      Read<T>(uncompressed.write()), is_little_endian);
}

void write_tag(std::ostream& stream, TagBase const* tag, Int space_dim) {
  if (is<I8>(tag)) {
    write_array(stream, tag->name(), tag->ncomps(), as<I8>(tag)->array());
  } else if (is<I32>(tag)) {
    write_array(stream, tag->name(), tag->ncomps(), as<I32>(tag)->array());
  } else if (is<I64>(tag)) {
    write_array(stream, tag->name(), tag->ncomps(), as<I64>(tag)->array());
  } else if (is<Real>(tag)) {
    Reals array = as<Real>(tag)->array();
    if (1 < space_dim && space_dim < 3) {
      if (tag->ncomps() == space_dim) {
        // VTK / ParaView expect vector fields to have 3 components
        // regardless of whether this is a 2D mesh or not.
        // this filter adds a 3rd zero component to any
        // fields with 2 components for 2D meshes
        write_array(
            stream, tag->name(), 3, resize_vectors(array, space_dim, 3));
      } else if (tag->ncomps() == symm_ncomps(space_dim)) {
        // Likewise, ParaView has component names specially set up for
        // 3D symmetric tensors
        write_array(stream, tag->name(), symm_ncomps(3),
            resize_symms(array, space_dim, 3));
      } else {
        write_array(stream, tag->name(), tag->ncomps(), array);
      }
    } else {
      write_array(stream, tag->name(), tag->ncomps(), array);
    }
  } else {
    Omega_h_fail("unknown tag type in write_tag");
  }
}

bool read_tag(std::istream& stream, Mesh* mesh, Int ent_dim,
    bool is_little_endian, bool is_compressed) {
  Omega_h_Type type = OMEGA_H_I8;
  std::string name;
  Int ncomps = -1;
  if (!read_array_start_tag(stream, &type, &name, &ncomps)) {
    return false;
  }
  /* tags like "global" are set by the construction mechanism,
     and it is somewhat complex to anticipate when they exist
     so we can just remove them if they are going to be reset. */
  mesh->remove_tag(ent_dim, name);
  auto size = mesh->nents(ent_dim) * ncomps;
  if (type == OMEGA_H_I8) {
    auto array = read_array<I8>(stream, size, is_little_endian, is_compressed);
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  } else if (type == OMEGA_H_I32) {
    auto array = read_array<I32>(stream, size, is_little_endian, is_compressed);
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  } else if (type == OMEGA_H_I64) {
    auto array = read_array<I64>(stream, size, is_little_endian, is_compressed);
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  } else {
    auto array =
        read_array<Real>(stream, size, is_little_endian, is_compressed);
    // undo the resizes done in write_tag()
    if (1 < mesh->dim() && mesh->dim() < 3) {
      if (ncomps == 3) {
        array = resize_vectors(array, 3, mesh->dim());
        ncomps = mesh->dim();
      } else if (ncomps == symm_ncomps(3)) {
        array = resize_symms(array, 3, mesh->dim());
        ncomps = symm_ncomps(mesh->dim());
      }
    }
    mesh->add_tag(ent_dim, name, ncomps, array, true);
  }
  auto et = xml::read_tag(stream);
  OMEGA_H_CHECK(et.elem_name == "DataArray");
  OMEGA_H_CHECK(et.type == xml::Tag::END);
  return true;
}

template <typename T>
Read<T> read_known_array(std::istream& stream, std::string const& name,
    LO nents, Int ncomps, bool is_little_endian, bool is_compressed) {
  auto st = xml::read_tag(stream);
  OMEGA_H_CHECK(st.elem_name == "DataArray");
  OMEGA_H_CHECK(st.type == xml::Tag::START);
  OMEGA_H_CHECK(st.attribs["Name"] == name);
  OMEGA_H_CHECK(st.attribs["type"] == Traits<T>::name());
  OMEGA_H_CHECK(st.attribs["NumberOfComponents"] == Omega_h::to_string(ncomps));
  auto array =
      read_array<T>(stream, nents * ncomps, is_little_endian, is_compressed);
  auto et = xml::read_tag(stream);
  OMEGA_H_CHECK(et.elem_name == "DataArray");
  OMEGA_H_CHECK(et.type == xml::Tag::END);
  return array;
}

enum {
  VTK_VERTEX = 1,
  VTK_POLY_VERTEX = 2,
  VTK_LINE = 3,
  VTK_POLY_LINE = 4,
  VTK_TRIANGLE = 5,
  VTK_TRIANGLE_STRIP = 6,
  VTK_POLYGON = 7,
  VTK_PIXEL = 8,
  VTK_QUAD = 9,
  VTK_TETRA = 10,
  VTK_VOXEL = 11,
  VTK_HEXAHEDRON = 12,
  VTK_WEDGE = 13,
  VTK_PYRAMID = 14
};

static I8 const vtk_types[DIMS] = {
    VTK_VERTEX, VTK_LINE, VTK_TRIANGLE, VTK_TETRA};

static void write_vtkfile_vtu_start_tag(std::ostream& stream) {
  stream << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"";
  if (is_little_endian_cpu())
    stream << "LittleEndian";
  else
    stream << "BigEndian";
  stream << "\" header_type=\"";
  stream << Traits<std::uint64_t>::name();
  stream << "\"";
#ifdef OMEGA_H_USE_ZLIB
  stream << " compressor=\"vtkZLibDataCompressor\"";
#endif
  stream << ">\n";
}

static void read_vtkfile_vtu_start_tag(
    std::istream& stream, bool* is_little_endian_out, bool* is_compressed_out) {
  auto st = xml::read_tag(stream);
  OMEGA_H_CHECK(st.elem_name == "VTKFile");
  OMEGA_H_CHECK(st.attribs["header_type"] == Traits<std::uint64_t>::name());
  auto is_little_endian = (st.attribs["byte_order"] == "LittleEndian");
  *is_little_endian_out = is_little_endian;
  auto is_compressed = (st.attribs.count("compressor") == 1);
  *is_compressed_out = is_compressed;
}

void write_piece_start_tag(
    std::ostream& stream, Mesh const* mesh, Int cell_dim) {
  stream << "<Piece NumberOfPoints=\"" << mesh->nverts() << "\"";
  stream << " NumberOfCells=\"" << mesh->nents(cell_dim) << "\">\n";
}

void read_piece_start_tag(
    std::istream& stream, LO* nverts_out, LO* ncells_out) {
  auto st = xml::read_tag(stream);
  OMEGA_H_CHECK(st.elem_name == "Piece");
  *nverts_out = std::stoi(st.attribs["NumberOfPoints"]);
  *ncells_out = std::stoi(st.attribs["NumberOfCells"]);
}

void write_connectivity(std::ostream& stream, Mesh* mesh, Int cell_dim) {
  Read<I8> types(mesh->nents(cell_dim), vtk_types[cell_dim]);
  write_array(stream, "types", 1, types);
  LOs ev2v = mesh->ask_verts_of(cell_dim);
  LOs ends(mesh->nents(cell_dim), simplex_degrees[cell_dim][VERT],
      simplex_degrees[cell_dim][VERT]);
  write_array(stream, "connectivity", 1, ev2v);
  write_array(stream, "offsets", 1, ends);
}

void read_connectivity(std::istream& stream, CommPtr comm, LO ncells,
    bool is_little_endian, bool is_compressed, Int* dim_out, LOs* ev2v_out) {
  auto types = read_known_array<I8>(
      stream, "types", ncells, 1, is_little_endian, is_compressed);
  Int dim = -1;
  if (types.size()) {
    auto type = types.get(0);
    if (type == VTK_TRIANGLE) dim = 2;
    if (type == VTK_TETRA) dim = 3;
  }
  dim = comm->allreduce(dim, OMEGA_H_MAX);
  OMEGA_H_CHECK(dim == 2 || dim == 3);
  *dim_out = dim;
  auto ev2v = read_known_array<LO>(stream, "connectivity", ncells * (dim + 1),
      1, is_little_endian, is_compressed);
  *ev2v_out = ev2v;
  read_known_array<LO>(
      stream, "offsets", ncells, 1, is_little_endian, is_compressed);
}

void write_locals(std::ostream& stream, Mesh* mesh, Int ent_dim) {
  write_array(stream, "local", 1, Read<LO>(mesh->nents(ent_dim), 0, 1));
}

void write_owners(std::ostream& stream, Mesh* mesh, Int ent_dim) {
  if (mesh->comm()->size() == 1) return;
  write_array(stream, "owner", 1, mesh->ask_owners(ent_dim).ranks);
}

void write_locals_and_owners(
    std::ostream& stream, Mesh* mesh, Int ent_dim, TagSet const& tags) {
  if (tags[size_t(ent_dim)].count("local")) {
    write_locals(stream, mesh, ent_dim);
  }
  if (tags[size_t(ent_dim)].count("owner")) {
    write_owners(stream, mesh, ent_dim);
  }
}

template <typename T>
void write_p_data_array(
    std::ostream& stream, std::string const& name, Int ncomps) {
  stream << "<PDataArray ";
  describe_array<T>(stream, name, ncomps);
  stream << "/>\n";
}

void write_p_data_array2(std::ostream& stream, std::string const& name,
    Int ncomps, Int Omega_h_Type) {
  switch (Omega_h_Type) {
    case OMEGA_H_I8:
      write_p_data_array<I8>(stream, name, ncomps);
      break;
    case OMEGA_H_I32:
      write_p_data_array<I32>(stream, name, ncomps);
      break;
    case OMEGA_H_I64:
      write_p_data_array<I64>(stream, name, ncomps);
      break;
    case OMEGA_H_F64:
      write_p_data_array<Real>(stream, name, ncomps);
      break;
  }
}

void write_p_tag(std::ostream& stream, TagBase const* tag, Int space_dim) {
  if (tag->type() == OMEGA_H_REAL) {
    if (1 < space_dim && space_dim < 3) {
      if (tag->ncomps() == space_dim) {
        write_p_data_array2(stream, tag->name(), 3, OMEGA_H_REAL);
      } else if (tag->ncomps() == symm_ncomps(space_dim)) {
        write_p_data_array2(stream, tag->name(), symm_ncomps(3), OMEGA_H_REAL);
      } else {
        write_p_data_array2(stream, tag->name(), tag->ncomps(), OMEGA_H_REAL);
      }
    } else {
      write_p_data_array2(stream, tag->name(), tag->ncomps(), OMEGA_H_REAL);
    }
  } else {
    write_p_data_array2(stream, tag->name(), tag->ncomps(), tag->type());
  }
}

std::string piece_filename(std::string const& piecepath, I32 rank) {
  return piecepath + '_' + Omega_h::to_string(rank) + ".vtu";
}

std::string get_rel_step_path(Int step) {
  return "steps/step_" + Omega_h::to_string(step);
}

std::string get_step_path(std::string const& root_path, Int step) {
  return root_path + '/' + get_rel_step_path(step);
}

}  // end anonymous namespace

std::string get_pvtu_path(std::string const& step_path) {
  return step_path + "/pieces.pvtu";
}

std::string get_pvd_path(std::string const& root_path) {
  return root_path + "/steps.pvd";
}

static void default_dim(Mesh* mesh, Int* cell_dim) {
  if (*cell_dim == -1) *cell_dim = mesh->dim();
}

void write_vtu(
    std::ostream& stream, Mesh* mesh, Int cell_dim, TagSet const& tags) {
  default_dim(mesh, &cell_dim);
  write_vtkfile_vtu_start_tag(stream);
  stream << "<UnstructuredGrid>\n";
  write_piece_start_tag(stream, mesh, cell_dim);
  stream << "<Cells>\n";
  write_connectivity(stream, mesh, cell_dim);
  stream << "</Cells>\n";
  stream << "<Points>\n";
  auto coords = mesh->coords();
  write_array(stream, "coordinates", 3, resize_vectors(coords, mesh->dim(), 3));
  stream << "</Points>\n";
  stream << "<PointData>\n";
  /* globals go first so read_vtu() knows where to find them */
  if (mesh->has_tag(VERT, "global") && tags[VERT].count("global")) {
    write_tag(stream, mesh->get_tag<GO>(VERT, "global"), mesh->dim());
  }
  write_locals_and_owners(stream, mesh, VERT, tags);
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tag = mesh->get_tag(VERT, i);
    if (tag->name() != "coordinates" && tag->name() != "global" &&
        tags[VERT].count(tag->name())) {
      write_tag(stream, tag, mesh->dim());
    }
  }
  stream << "</PointData>\n";
  stream << "<CellData>\n";
  /* globals go first so read_vtu() knows where to find them */
  if (mesh->has_tag(cell_dim, "global") &&
      tags[size_t(cell_dim)].count("global")) {
    write_tag(stream, mesh->get_tag<GO>(cell_dim, "global"), mesh->dim());
  }
  write_locals_and_owners(stream, mesh, cell_dim, tags);
  for (Int i = 0; i < mesh->ntags(cell_dim); ++i) {
    auto tag = mesh->get_tag(cell_dim, i);
    if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
      write_tag(stream, tag, mesh->dim());
    }
  }
  stream << "</CellData>\n";
  stream << "</Piece>\n";
  stream << "</UnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

void read_vtu(std::istream& stream, CommPtr comm, Mesh* mesh) {
  mesh->set_comm(comm);
  mesh->set_parting(OMEGA_H_ELEM_BASED);
  read_vtu_ents(stream, mesh);
}

void read_vtu_ents(std::istream& stream, Mesh* mesh) {
  bool is_little_endian, is_compressed;
  read_vtkfile_vtu_start_tag(stream, &is_little_endian, &is_compressed);
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "UnstructuredGrid");
  LO nverts, ncells;
  read_piece_start_tag(stream, &nverts, &ncells);
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "Cells");
  auto comm = mesh->comm();
  Int dim;
  LOs ev2v;
  read_connectivity(
      stream, comm, ncells, is_little_endian, is_compressed, &dim, &ev2v);
  mesh->set_dim(dim);
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "Cells");
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "Points");
  auto coords = read_known_array<Real>(
      stream, "coordinates", nverts, 3, is_little_endian, is_compressed);
  if (dim < 3) coords = resize_vectors(coords, 3, dim);
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "Points");
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "PointData");
  GOs vert_globals;
  if (mesh->could_be_shared(VERT)) {
    vert_globals = read_known_array<GO>(
        stream, "global", nverts, 1, is_little_endian, is_compressed);
  } else {
    vert_globals = Read<GO>(nverts, 0, 1);
  }
  build_verts_from_globals(mesh, vert_globals);
  mesh->add_tag(VERT, "coordinates", dim, coords, true);
  while (read_tag(stream, mesh, VERT, is_little_endian, is_compressed))
    ;
  mesh->remove_tag(VERT, "local");
  mesh->remove_tag(VERT, "owner");
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "CellData");
  GOs elem_globals;
  if (mesh->could_be_shared(dim)) {
    elem_globals = read_known_array<GO>(
        stream, "global", ncells, 1, is_little_endian, is_compressed);
  } else {
    elem_globals = Read<GO>(ncells, 0, 1);
  }
  build_ents_from_elems2verts(mesh, ev2v, vert_globals, elem_globals);
  while (read_tag(stream, mesh, dim, is_little_endian, is_compressed))
    ;
  mesh->remove_tag(dim, "local");
  mesh->remove_tag(dim, "owner");
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "Piece");
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "UnstructuredGrid");
  OMEGA_H_CHECK(xml::read_tag(stream).elem_name == "VTKFile");
}

void write_vtu(
    std::string const& filename, Mesh* mesh, Int cell_dim, TagSet const& tags) {
  std::ofstream file(filename.c_str());
  OMEGA_H_CHECK(file.is_open());
  ask_for_mesh_tags(mesh, tags);
  write_vtu(file, mesh, cell_dim, tags);
}

void write_vtu(std::string const& filename, Mesh* mesh, Int cell_dim) {
  write_vtu(filename, mesh, cell_dim, get_all_vtk_tags(mesh));
}

void write_vtu(std::string const& filename, Mesh* mesh) {
  write_vtu(filename, mesh, mesh->dim());
}

void write_pvtu(std::ostream& stream, Mesh* mesh, Int cell_dim,
    std::string const& piecepath, TagSet const& tags) {
  ask_for_mesh_tags(mesh, tags);
  stream << "<VTKFile type=\"PUnstructuredGrid\">\n";
  stream << "<PUnstructuredGrid";
  stream << " GhostLevel=\"";
  if (mesh->parting() == OMEGA_H_GHOSTED) {
    stream << mesh->nghost_layers();
  } else {
    if (mesh->parting() == OMEGA_H_VERT_BASED && can_print(mesh) == 0) {
      std::cerr << "WARNING: a pvtu file may be written from a "
                   "vertex-partitioned mesh, but NOT read back in\n";
    }
    stream << Int(0);
  }
  stream << "\">\n";
  stream << "<PPoints>\n";
  write_p_data_array<Real>(stream, "coordinates", 3);
  stream << "</PPoints>\n";
  stream << "<PPointData>\n";
  if (mesh->has_tag(VERT, "global") && tags[VERT].count("global")) {
    write_p_tag(stream, mesh->get_tag<GO>(VERT, "global"), mesh->dim());
  }
  if (tags[0].count("local")) {
    write_p_data_array2(stream, "local", 1, OMEGA_H_I32);
  }
  if (mesh->comm()->size() > 1 && tags[0].count("owner")) {
    write_p_data_array2(stream, "owner", 1, OMEGA_H_I32);
  }
  for (Int i = 0; i < mesh->ntags(VERT); ++i) {
    auto tag = mesh->get_tag(VERT, i);
    if (tag->name() != "coordinates" && tag->name() != "global" &&
        tags[VERT].count(tag->name())) {
      write_p_tag(stream, tag, mesh->dim());
    }
  }
  stream << "</PPointData>\n";
  stream << "<PCellData>\n";
  if (mesh->has_tag(cell_dim, "global") && tags[size_t(cell_dim)].count("global")) {
    write_p_tag(stream, mesh->get_tag<GO>(cell_dim, "global"), mesh->dim());
  }
  if (tags[size_t(cell_dim)].count("local")) {
    write_p_data_array2(stream, "local", 1, OMEGA_H_I32);
  }
  if (mesh->comm()->size() > 1 && tags[size_t(cell_dim)].count("owner")) {
    write_p_data_array2(stream, "owner", 1, OMEGA_H_I32);
  }
  for (Int i = 0; i < mesh->ntags(cell_dim); ++i) {
    auto tag = mesh->get_tag(cell_dim, i);
    if (tag->name() != "global" && tags[size_t(cell_dim)].count(tag->name())) {
      write_p_tag(stream, tag, mesh->dim());
    }
  }
  stream << "</PCellData>\n";
  for (I32 i = 0; i < mesh->comm()->size(); ++i) {
    stream << "<Piece Source=\"" << piece_filename(piecepath, i) << "\"/>\n";
  }
  stream << "</PUnstructuredGrid>\n";
  stream << "</VTKFile>\n";
}

void write_pvtu(std::string const& filename, Mesh* mesh, Int cell_dim,
    std::string const& piecepath, TagSet const& tags) {
  std::ofstream file(filename.c_str());
  OMEGA_H_CHECK(file.is_open());
  write_pvtu(file, mesh, cell_dim, piecepath, tags);
}

void read_pvtu(std::istream& stream, CommPtr comm, I32* npieces_out,
    std::string* vtupath_out, Int* nghost_layers_out) {
  I32 npieces = 0;
  std::string vtupath;
  Int nghost_layers = 0;
  for (std::string line; std::getline(stream, line);) {
    xml::Tag tag;
    if (!xml::parse_tag(line, &tag)) continue;
    if (tag.elem_name == "Piece") {
      if (npieces == comm->rank()) {
        vtupath = tag.attribs["Source"];
      }
      ++npieces;
    }
    if (tag.elem_name == "PUnstructuredGrid") {
      if (tag.attribs.find("GhostLevel") != tag.attribs.end()) {
        auto& ghostlevelstr = tag.attribs["GhostLevel"];
        nghost_layers = std::atoi(ghostlevelstr.c_str());
      }
    }
  }
  OMEGA_H_CHECK(npieces >= 1);
  OMEGA_H_CHECK(npieces <= comm->size());
  *npieces_out = npieces;
  *vtupath_out = vtupath;
  *nghost_layers_out = nghost_layers;
}

void read_pvtu(std::string const& pvtupath, CommPtr comm, I32* npieces_out,
    std::string* vtupath_out, Int* nghost_layers_out) {
  auto parentpath = parent_path(pvtupath);
  std::string vtupath;
  std::ifstream stream(pvtupath.c_str());
  if (!stream.is_open()) {
    Omega_h_fail("couldn't open \"%s\"\n", pvtupath.c_str());
  }
  read_pvtu(stream, comm, npieces_out, &vtupath, nghost_layers_out);
  vtupath = parentpath + "/" + vtupath;
  *vtupath_out = vtupath;
}

void write_parallel(
    std::string const& path, Mesh* mesh, Int cell_dim, TagSet const& tags) {
  default_dim(mesh, &cell_dim);
  ask_for_mesh_tags(mesh, tags);
  auto rank = mesh->comm()->rank();
  if (rank == 0) {
    safe_mkdir(path.c_str());
  }
  mesh->comm()->barrier();
  auto piecesdir = path + "/pieces";
  if (rank == 0) {
    safe_mkdir(piecesdir.c_str());
  }
  mesh->comm()->barrier();
  auto piecepath = piecesdir + "/piece";
  auto pvtuname = get_pvtu_path(path);
  if (rank == 0) {
    write_pvtu(pvtuname, mesh, cell_dim, "pieces/piece", tags);
  }
  write_vtu(piece_filename(piecepath, rank), mesh, cell_dim, tags);
}

void write_parallel(std::string const& path, Mesh* mesh, Int cell_dim) {
  write_parallel(path, mesh, cell_dim, get_all_vtk_tags(mesh));
}

void read_parallel(std::string const& pvtupath, CommPtr comm, Mesh* mesh) {
  I32 npieces;
  std::string vtupath;
  Int nghost_layers;
  read_pvtu(pvtupath, comm, &npieces, &vtupath, &nghost_layers);
  bool in_subcomm = (comm->rank() < npieces);
  auto subcomm = comm->split(I32(!in_subcomm), 0);
  if (in_subcomm) {
    std::ifstream vtustream(vtupath.c_str());
    OMEGA_H_CHECK(vtustream.is_open());
    mesh->set_comm(subcomm);
    if (nghost_layers == 0) {
      mesh->set_parting(OMEGA_H_ELEM_BASED, 0, false);
    } else {
      mesh->set_parting(OMEGA_H_GHOSTED, nghost_layers, false);
    }
    read_vtu_ents(vtustream, mesh);
  }
  mesh->set_comm(comm);
}

std::streampos write_initial_pvd(std::string const& root_path) {
  std::string pvdpath = get_pvd_path(root_path);
  std::ofstream file(pvdpath.c_str());
  OMEGA_H_CHECK(file.is_open());
  file << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
  file << "<Collection>\n";
  auto pos = file.tellp();
  file << "</Collection>\n";
  file << "</VTKFile>\n";
  return pos;
}

void update_pvd(std::string const& root_path, std::streampos* pos_inout,
    Int step, Real time) {
  std::string pvdpath = get_pvd_path(root_path);
  std::fstream file;
  file.open(pvdpath.c_str(), std::ios::out | std::ios::in);
  OMEGA_H_CHECK(file.is_open());
  file.seekp(*pos_inout);
  file << "<DataSet timestep=\"" << time << "\" part=\"0\" ";
  auto relstep = get_rel_step_path(step);
  auto relpvtu = get_pvtu_path(relstep);
  file << "file=\"" << relpvtu << "\"/>\n";
  *pos_inout = file.tellp();
  file << "</Collection>\n";
  file << "</VTKFile>\n";
}

void read_pvd(std::istream& stream, std::vector<Real>* times_out,
    std::vector<std::string>* pvtupaths_out) {
  std::vector<Real> times;
  std::vector<std::string> pvtupaths;
  for (std::string line; std::getline(stream, line);) {
    xml::Tag tag;
    if (!xml::parse_tag(line, &tag)) continue;
    if (tag.elem_name != "DataSet") continue;
    times.push_back(std::stod(tag.attribs["timestep"]));
    pvtupaths.push_back(tag.attribs["file"]);
  }
  *times_out = times;
  *pvtupaths_out = pvtupaths;
}

void read_pvd(std::string const& pvdpath, std::vector<Real>* times_out,
    std::vector<std::string>* pvtupaths_out) {
  std::vector<Real> times;
  std::vector<std::string> pvtupaths;
  std::ifstream pvdstream(pvdpath.c_str());
  OMEGA_H_CHECK(pvdstream.is_open());
  read_pvd(pvdstream, &times, &pvtupaths);
  auto parentpath = parent_path(pvdpath);
  for (auto& pvtupath : pvtupaths) pvtupath = parentpath + "/" + pvtupath;
  *times_out = times;
  *pvtupaths_out = pvtupaths;
}

Writer::Writer()
    : mesh_(nullptr),
      root_path_("/not-set"),
      cell_dim_(-1),
      step_(-1),
      pvd_pos_(0) {}

Writer::Writer(Writer const& other)
    : mesh_(other.mesh_),
      root_path_(other.root_path_),
      cell_dim_(other.cell_dim_),
      step_(other.step_),
      pvd_pos_(other.pvd_pos_) {}

Writer& Writer::operator=(Writer const& other) {
  mesh_ = other.mesh_;
  root_path_ = other.root_path_;
  cell_dim_ = other.cell_dim_;
  step_ = other.step_;
  pvd_pos_ = other.pvd_pos_;
  return *this;
}

Writer::~Writer() {}

Writer::Writer(std::string const& root_path, Mesh* mesh, Int cell_dim)
    : mesh_(mesh),
      root_path_(root_path),
      cell_dim_(cell_dim),
      step_(0),
      pvd_pos_(0) {
  default_dim(mesh_, &cell_dim_);
  auto comm = mesh->comm();
  auto rank = comm->rank();
  if (rank == 0) safe_mkdir(root_path_.c_str());
  comm->barrier();
  auto stepsdir = root_path_ + "/steps";
  if (rank == 0) safe_mkdir(stepsdir.c_str());
  comm->barrier();
  if (rank == 0) {
    pvd_pos_ = write_initial_pvd(root_path);
  }
}

void Writer::write(Real time, TagSet const& tags) {
  write_parallel(get_step_path(root_path_, step_), mesh_, cell_dim_, tags);
  if (mesh_->comm()->rank() == 0) {
    update_pvd(root_path_, &pvd_pos_, step_, time);
  }
  ++step_;
}

void Writer::write(Real time) { this->write(time, get_all_vtk_tags(mesh_)); }

void Writer::write() { this->write(Real(step_)); }

FullWriter::FullWriter() {}

FullWriter::FullWriter(FullWriter const& other) : writers_(other.writers_) {}

FullWriter& FullWriter::operator=(FullWriter const& other) {
  writers_ = other.writers_;
  return *this;
}

FullWriter::~FullWriter() {}

FullWriter::FullWriter(std::string const& root_path, Mesh* mesh) {
  auto comm = mesh->comm();
  auto rank = comm->rank();
  if (rank == 0) safe_mkdir(root_path.c_str());
  comm->barrier();
  for (Int i = EDGE; i <= mesh->dim(); ++i)
    writers_.push_back(Writer(root_path + "/" + plural_names[i], mesh, i));
}

void FullWriter::write(Real time) {
  for (auto& writer : writers_) writer.write(time);
}

void FullWriter::write() {
  for (auto& writer : writers_) writer.write();
}

}  // end namespace vtk

}  // end namespace Omega_h
