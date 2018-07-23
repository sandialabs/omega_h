#include "Omega_h_file.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>

#ifdef OMEGA_H_USE_ZLIB
#include <zlib.h>
#endif

#include "Omega_h_array_ops.hpp"
#include "Omega_h_inertia.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

bool ends_with(std::string const& s, std::string const& suffix) {
  if (s.length() < suffix.length()) return false;
  return 0 == s.compare(s.length() - suffix.length(), suffix.length(), suffix);
}

bool is_little_endian_cpu() {
  static std::uint16_t const endian_canary = 0x1;
  std::uint8_t const* p = reinterpret_cast<std::uint8_t const*>(&endian_canary);
  return *p == 0x1;
}

void safe_mkdir(const char* path) {
  OMEGA_H_TIME_FUNCTION;
  mode_t const mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST) {
    Omega_h_fail(
        "omega_h could not create directory \"%s\", got error \"%s\"\n", path,
        std::strerror(errno));
  }
}

bool file_exists(const char* path) {
  struct stat info;
  if (stat(path, &info) != 0) return false;
  OMEGA_H_CHECK(info.st_mode & S_IFREG);
  return true;
}

bool directory_exists(const char* path) {
  struct stat info;
  if (stat(path, &info) != 0) return false;
  OMEGA_H_CHECK(info.st_mode & S_IFDIR);
  return true;
}

std::string parent_path(std::string const& path) {
  auto pos = path.find_last_of('/');
  OMEGA_H_CHECK(pos != std::string::npos);
  return path.substr(0, pos);
}

std::string path_leaf_name(std::string const& path) {
  auto pos = path.find_last_of('/');
  if (pos == std::string::npos) return path;
  return path.substr(pos + 1, std::string::npos);
}

namespace binary {

namespace {

static_assert(sizeof(Int) == 4, "osh format assumes 32 bit Int");
static_assert(sizeof(LO) == 4, "osh format assumes 32 bit LO");
static_assert(sizeof(GO) == 8, "osh format assumes 64 bit GO");
static_assert(sizeof(Real) == 8, "osh format assumes 64 bit Real");

OMEGA_H_INLINE std::uint32_t bswap32(std::uint32_t a) {
#ifdef OMEGA_H_USE_CUDA
  a = ((a & 0x000000FF) << 24) | ((a & 0x0000FF00) << 8) |
      ((a & 0x00FF0000) >> 8) | ((a & 0xFF000000) >> 24);
#else
  a = __builtin_bswap32(a);
#endif
  return a;
}

OMEGA_H_INLINE std::uint64_t bswap64(std::uint64_t a) {
#ifdef OMEGA_H_USE_CUDA
  a = ((a & 0x00000000000000FFULL) << 56) |
      ((a & 0x000000000000FF00ULL) << 40) |
      ((a & 0x0000000000FF0000ULL) << 24) | ((a & 0x00000000FF000000ULL) << 8) |
      ((a & 0x000000FF00000000ULL) >> 8) | ((a & 0x0000FF0000000000ULL) >> 24) |
      ((a & 0x00FF000000000000ULL) >> 40) | ((a & 0xFF00000000000000ULL) >> 56);
#else
  a = __builtin_bswap64(a);
#endif
  return a;
}

template <typename T, size_t size = sizeof(T)>
struct SwapBytes;

template <typename T>
struct SwapBytes<T, 1> {
  OMEGA_H_INLINE static void swap(T*) {}
};

template <typename T>
struct SwapBytes<T, 4> {
  OMEGA_H_INLINE static void swap(T* ptr) {
    std::uint32_t* p2 = reinterpret_cast<std::uint32_t*>(ptr);
    *p2 = bswap32(*p2);
  }
};

template <typename T>
struct SwapBytes<T, 8> {
  OMEGA_H_INLINE static void swap(T* ptr) {
    std::uint64_t* p2 = reinterpret_cast<std::uint64_t*>(ptr);
    *p2 = bswap64(*p2);
  }
};

template <typename T>
OMEGA_H_INLINE void swap_bytes(T* ptr) {
  SwapBytes<T>::swap(ptr);
}

unsigned char const magic[2] = {0xa1, 0x1a};

}  // end anonymous namespace

template <typename T>
void swap_if_needed(T& val, bool is_little_endian) {
  if (is_little_endian != is_little_endian_cpu()) {
    swap_bytes(&val);
  }
}

template <typename T>
Read<T> swap_if_needed(Read<T> array, bool is_little_endian) {
  if (is_little_endian == is_little_endian_cpu()) {
    return array;
  }
  Write<T> out = deep_copy(array);
  auto f = OMEGA_H_LAMBDA(LO i) { swap_bytes(&out[i]); };
  parallel_for(out.size(), f, "swap_if_needed");
  return out;
}

template <typename T>
void write_value(std::ostream& stream, T val) {
  swap_if_needed(val);
  stream.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
void read_value(std::istream& stream, T& val) {
  stream.read(reinterpret_cast<char*>(&val), sizeof(T));
  swap_if_needed(val);
}

template <typename T>
void write_array(std::ostream& stream, Read<T> array) {
  LO size = array.size();
  write_value(stream, size);
  Read<T> swapped = swap_if_needed(array, true);
  HostRead<T> uncompressed(swapped);
  I64 uncompressed_bytes =
      static_cast<I64>(static_cast<std::size_t>(size) * sizeof(T));
#ifdef OMEGA_H_USE_ZLIB
  uLong source_bytes = static_cast<uLong>(uncompressed_bytes);
  uLong dest_bytes = ::compressBound(source_bytes);
  auto compressed = new ::Bytef[dest_bytes];
  int ret = ::compress2(compressed, &dest_bytes,
      reinterpret_cast<const ::Bytef*>(nonnull(uncompressed.data())),
      source_bytes, Z_BEST_SPEED);
  OMEGA_H_CHECK(ret == Z_OK);
  I64 compressed_bytes = static_cast<I64>(dest_bytes);
  write_value(stream, compressed_bytes);
  stream.write(reinterpret_cast<const char*>(compressed), compressed_bytes);
  delete[] compressed;
#else
  stream.write(reinterpret_cast<const char*>(nonnull(uncompressed.data())),
      uncompressed_bytes);
#endif
}

template <typename T>
void read_array(std::istream& stream, Read<T>& array, bool is_compressed) {
  LO size;
  read_value(stream, size);
  OMEGA_H_CHECK(size >= 0);
  I64 uncompressed_bytes =
      static_cast<I64>(static_cast<std::size_t>(size) * sizeof(T));
  HostWrite<T> uncompressed(size);
#ifdef OMEGA_H_USE_ZLIB
  if (is_compressed) {
    I64 compressed_bytes;
    read_value(stream, compressed_bytes);
    OMEGA_H_CHECK(compressed_bytes >= 0);
    auto compressed = new ::Bytef[compressed_bytes];
    stream.read(reinterpret_cast<char*>(compressed), compressed_bytes);
    uLong dest_bytes = static_cast<uLong>(uncompressed_bytes);
    uLong source_bytes = static_cast<uLong>(compressed_bytes);
    ::Bytef* uncompressed_ptr =
        reinterpret_cast< ::Bytef*>(nonnull(uncompressed.data()));
    int ret =
        ::uncompress(uncompressed_ptr, &dest_bytes, compressed, source_bytes);
    OMEGA_H_CHECK(ret == Z_OK);
    OMEGA_H_CHECK(dest_bytes == static_cast<uLong>(uncompressed_bytes));
    delete[] compressed;
  } else
#else
  OMEGA_H_CHECK(is_compressed == false);
#endif
  {
    stream.read(reinterpret_cast<char*>(nonnull(uncompressed.data())),
        uncompressed_bytes);
  }
  array = swap_if_needed(Read<T>(uncompressed.write()), true);
}

void write(std::ostream& stream, std::string const& val) {
  I32 len = static_cast<I32>(val.length());
  write_value(stream, len);
  stream.write(val.c_str(), len);
}

void read(std::istream& stream, std::string& val) {
  I32 len;
  read_value(stream, len);
  OMEGA_H_CHECK(len >= 0);
  val.resize(static_cast<std::size_t>(len));
  stream.read(&val[0], len);
}

static void write_meta(std::ostream& stream, Mesh const* mesh) {
  auto family = I8(mesh->family());
  write_value(stream, family);
  auto dim = I8(mesh->dim());
  write_value(stream, dim);
  I32 comm_size = mesh->comm()->size();
  write_value(stream, comm_size);
  I32 comm_rank = mesh->comm()->rank();
  write_value(stream, comm_rank);
  I8 parting = mesh->parting();
  write_value(stream, parting);
  I32 nghost_layers = mesh->nghost_layers();
  write_value(stream, nghost_layers);
  auto hints = mesh->rib_hints();
  I8 have_hints = (hints != nullptr);
  write_value(stream, have_hints);
  if (have_hints) {
    auto naxes = I32(hints->axes.size());
    write_value(stream, naxes);
    for (auto axis : hints->axes) {
      for (Int i = 0; i < 3; ++i) write_value(stream, axis[i]);
    }
  }
}

static void read_meta(std::istream& stream, Mesh* mesh, Int version) {
  if (version >= 7) {
    I8 family;
    read_value(stream, family);
    mesh->set_family(Omega_h_Family(family));
  }
  I8 dim;
  read_value(stream, dim);
  mesh->set_dim(Int(dim));
  I32 comm_size;
  read_value(stream, comm_size);
  OMEGA_H_CHECK(mesh->comm()->size() == comm_size);
  I32 comm_rank;
  read_value(stream, comm_rank);
  OMEGA_H_CHECK(mesh->comm()->rank() == comm_rank);
  I8 parting_i8;
  read_value(stream, parting_i8);
  OMEGA_H_CHECK(parting_i8 == I8(OMEGA_H_ELEM_BASED) ||
                parting_i8 == I8(OMEGA_H_GHOSTED) ||
                parting_i8 == I8(OMEGA_H_VERT_BASED));
  if (version >= 3) {
    I32 nghost_layers;
    read_value(stream, nghost_layers);
    mesh->set_parting(Omega_h_Parting(parting_i8), nghost_layers, false);
  } else {
    mesh->set_parting(Omega_h_Parting(parting_i8));
  }
  I8 have_hints;
  read_value(stream, have_hints);
  if (have_hints) {
    I32 naxes;
    read_value(stream, naxes);
    auto hints = std::make_shared<inertia::Rib>();
    for (I32 i = 0; i < naxes; ++i) {
      Vector<3> axis;
      for (Int j = 0; j < 3; ++j) read_value(stream, axis[j]);
      hints->axes.push_back(axis);
    }
    mesh->set_rib_hints(hints);
  }
  if (version < 6) {
    I8 keeps_canon;
    read_value(stream, keeps_canon);
  }
}

static void write_tag(std::ostream& stream, TagBase const* tag) {
  std::string name = tag->name();
  write(stream, name);
  auto ncomps = I8(tag->ncomps());
  write_value(stream, ncomps);
  I8 type = tag->type();
  write_value(stream, type);
  if (is<I8>(tag)) {
    write_array(stream, as<I8>(tag)->array());
  } else if (is<I32>(tag)) {
    write_array(stream, as<I32>(tag)->array());
  } else if (is<I64>(tag)) {
    write_array(stream, as<I64>(tag)->array());
  } else if (is<Real>(tag)) {
    write_array(stream, as<Real>(tag)->array());
  } else {
    Omega_h_fail("unexpected tag type in binary write\n");
  }
}

static void read_tag(
    std::istream& stream, Mesh* mesh, Int d, bool is_compressed, I32 version) {
  std::string name;
  read(stream, name);
  I8 ncomps;
  read_value(stream, ncomps);
  I8 type;
  read_value(stream, type);
  if (version < 5) {
    I8 xfer_i8;
    read_value(stream, xfer_i8);
    if (2 <= version) {
      I8 outflags_i8;
      read_value(stream, outflags_i8);
    }
  }
  if (type == OMEGA_H_I8) {
    Read<I8> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, array, true);
  } else if (type == OMEGA_H_I32) {
    Read<I32> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, array, true);
  } else if (type == OMEGA_H_I64) {
    Read<I64> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, array, true);
  } else if (type == OMEGA_H_F64) {
    Read<Real> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, array, true);
  } else {
    Omega_h_fail("unexpected tag type in binary read\n");
  }
}

static void write_sets(std::ostream& stream, Mesh* mesh) {
  auto n = I32(mesh->class_sets.size());
  write_value(stream, n);
  for (auto& set : mesh->class_sets) {
    auto& name = set.first;
    write(stream, name);
    auto npairs = I32(set.second.size());
    write_value(stream, npairs);
    for (auto& pair : set.second) {
      write_value(stream, pair.dim);
      write_value(stream, pair.id);
    }
  }
}

static void read_sets(std::istream& stream, Mesh* mesh) {
  I32 n;
  read_value(stream, n);
  for (I32 i = 0; i < n; ++i) {
    std::string name;
    read(stream, name);
    I32 npairs;
    read_value(stream, npairs);
    for (I32 j = 0; j < npairs; ++j) {
      ClassPair pair;
      read_value(stream, pair.dim);
      read_value(stream, pair.id);
      mesh->class_sets[name].push_back(pair);
    }
  }
}

void write(std::ostream& stream, Mesh* mesh) {
  begin_code("binary::write(stream,Mesh)");
  stream.write(reinterpret_cast<const char*>(magic), sizeof(magic));
// write_value(stream, latest_version); moved to /version at version 4
#ifdef OMEGA_H_USE_ZLIB
  I8 is_compressed = true;
#else
  I8 is_compressed = false;
#endif
  write_value(stream, is_compressed);
  write_meta(stream, mesh);
  LO nverts = mesh->nverts();
  write_value(stream, nverts);
  for (Int d = 1; d <= mesh->dim(); ++d) {
    auto down = mesh->ask_down(d, d - 1);
    write_array(stream, down.ab2b);
    if (d > 1) {
      write_array(stream, down.codes);
    }
  }
  for (Int d = 0; d <= mesh->dim(); ++d) {
    auto nsaved_tags = mesh->ntags(d);
    write_value(stream, nsaved_tags);
    for (Int i = 0; i < mesh->ntags(d); ++i) {
      write_tag(stream, mesh->get_tag(d, i));
    }
    if (mesh->comm()->size() > 1) {
      auto owners = mesh->ask_owners(d);
      write_array(stream, owners.ranks);
      write_array(stream, owners.idxs);
    }
  }
  write_sets(stream, mesh);
  end_code();
}

void read(std::istream& stream, Mesh* mesh, I32 version) {
  unsigned char magic_in[2];
  stream.read(reinterpret_cast<char*>(magic_in), sizeof(magic));
  OMEGA_H_CHECK(magic_in[0] == magic[0]);
  OMEGA_H_CHECK(magic_in[1] == magic[1]);
  if (version == -1) read_value(stream, version);
  OMEGA_H_CHECK(version >= 1);
  OMEGA_H_CHECK(version <= latest_version);
  I8 is_compressed;
  read_value(stream, is_compressed);
#ifndef OMEGA_H_USE_ZLIB
  OMEGA_H_CHECK(!is_compressed);
#endif
  read_meta(stream, mesh, version);
  LO nverts;
  read_value(stream, nverts);
  mesh->set_verts(nverts);
  for (Int d = 1; d <= mesh->dim(); ++d) {
    Adj down;
    read_array(stream, down.ab2b, is_compressed);
    if (d > 1) {
      read_array(stream, down.codes, is_compressed);
    }
    mesh->set_ents(d, down);
  }
  for (Int d = 0; d <= mesh->dim(); ++d) {
    Int ntags;
    read_value(stream, ntags);
    for (Int i = 0; i < ntags; ++i) {
      read_tag(stream, mesh, d, is_compressed, version);
    }
    if (mesh->comm()->size() > 1) {
      Remotes owners;
      read_array(stream, owners.ranks, is_compressed);
      read_array(stream, owners.idxs, is_compressed);
      mesh->set_owners(d, owners);
    }
  }
  if (version >= 8) {
    read_sets(stream, mesh);
  }
}

static void write_int_file(std::string const& filepath, Mesh* mesh, I32 value) {
  if (mesh->comm()->rank() == 0) {
    std::ofstream file(filepath.c_str());
    OMEGA_H_CHECK(file.is_open());
    file << value << '\n';
  }
}

static void write_nparts(std::string const& path, Mesh* mesh) {
  write_int_file(path + "/nparts", mesh, mesh->comm()->size());
}

static void write_version(std::string const& path, Mesh* mesh) {
  write_int_file(path + "/version", mesh, latest_version);
}

I32 read_nparts(std::string const& path, CommPtr comm) {
  I32 nparts;
  if (comm->rank() == 0) {
    auto filepath = path + "/nparts";
    std::ifstream file(filepath.c_str());
    if (!file.is_open()) {
      Omega_h_fail("could not open file \"%s\"\n", filepath.c_str());
    }
    file >> nparts;
    if (!file) {
      Omega_h_fail("could not read file \"%s\"\n", filepath.c_str());
    }
  }
  comm->bcast(nparts);
  return nparts;
}

I32 read_version(std::string const& path, CommPtr comm) {
  I32 version;
  if (comm->rank() == 0) {
    auto filepath = path + "/version";
    std::ifstream file(filepath.c_str());
    if (!file.is_open()) {
      version = -1;
    } else {
      file >> version;
      if (!file) {
        Omega_h_fail("could not read file \"%s\"\n", filepath.c_str());
      }
    }
  }
  comm->bcast(version);
  return version;
}

void write(std::string const& path, Mesh* mesh) {
  begin_code("binary::write(path,Mesh)");
  if (!ends_with(path, ".osh") && can_print(mesh)) {
    std::cout
        << "it is strongly recommended to end Omega_h paths in \".osh\",\n";
    std::cout << "instead of just \"" << path << "\"\n";
  }
  safe_mkdir(path.c_str());
  mesh->comm()->barrier();
  auto filepath = path + "/" + to_string(mesh->comm()->rank()) + ".osh";
  std::ofstream file(filepath.c_str());
  OMEGA_H_CHECK(file.is_open());
  write(file, mesh);
  write_nparts(path, mesh);
  write_version(path, mesh);
  mesh->comm()->barrier();
  end_code();
}

void read_in_comm(
    std::string const& path, CommPtr comm, Mesh* mesh, I32 version) {
  mesh->set_comm(comm);
  auto filepath = path + "/" + to_string(mesh->comm()->rank());
  if (version != -1) filepath += ".osh";
  std::ifstream file(filepath.c_str());
  OMEGA_H_CHECK(file.is_open());
  read(file, mesh, version);
}

I32 read(std::string const& path, CommPtr comm, Mesh* mesh, bool strict) {
  auto nparts = read_nparts(path, comm);
  auto version = read_version(path, comm);
  if (strict) {
    if (nparts != comm->size()) {
      Omega_h_fail(
          "Mesh \"%s\" is being read in strict mode"
          " (no repartitioning) and its number of parts %d"
          " doesn't match the number of MPI ranks %d\n",
          path.c_str(), nparts, comm->size());
    }
    read_in_comm(path, comm, mesh, version);
  } else {
    if (nparts > comm->size()) {
      Omega_h_fail(
          "path \"%s\" contains %d parts, but only %d ranks are reading it\n",
          path.c_str(), nparts, comm->size());
    }
    auto in_subcomm = (comm->rank() < nparts);
    auto subcomm = comm->split(I32(!in_subcomm), 0);
    if (in_subcomm) {
      read_in_comm(path, subcomm, mesh, version);
    }
    mesh->set_comm(comm);
  }
  return nparts;
}

Mesh read(std::string const& path, Library* lib, bool strict) {
  return binary::read(path, lib->world(), strict);
}

Mesh read(std::string const& path, CommPtr comm, bool strict) {
  auto mesh = Mesh(comm->library());
  binary::read(path, comm, &mesh, strict);
  return mesh;
}

#define OMEGA_H_INST(T)                                                        \
  template void swap_if_needed(T& val, bool is_little_endian);                 \
  template Read<T> swap_if_needed(Read<T> array, bool is_little_endian);       \
  template void write_value(std::ostream& stream, T val);                      \
  template void read_value(std::istream& stream, T& val);                      \
  template void write_array(std::ostream& stream, Read<T> array);              \
  template void read_array(                                                    \
      std::istream& stream, Read<T>& array, bool is_compressed);
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

// for VTK compression headers
template void swap_if_needed(std::uint64_t& val, bool is_little_endian);

}  // end namespace binary

void write_reals_txt(std::string const& filename, Reals a, Int ncomps) {
  std::ofstream stream(filename.c_str());
  write_reals_txt(stream, a, ncomps);
}

void write_reals_txt(std::ostream& stream, Reals a, Int ncomps) {
  auto n = divide_no_remainder(a.size(), ncomps);
  auto h_a = HostRead<Real>(a);
  stream << std::scientific << std::setprecision(17);
  for (LO i = 0; i < n; ++i) {
    for (Int j = 0; j < ncomps; ++j) {
      stream << h_a[i * ncomps + j];
      if (j < ncomps - 1) stream << ' ';
    }
    stream << '\n';
  }
}

Reals read_reals_txt(std::string const& filename, LO n, Int ncomps) {
  std::ifstream stream(filename.c_str());
  return read_reals_txt(stream, n, ncomps);
}

Reals read_reals_txt(std::istream& stream, LO n, Int ncomps) {
  auto h_a = HostWrite<Real>(n * ncomps);
  for (LO i = 0; i < n; ++i) {
    for (Int j = 0; j < ncomps; ++j) {
      stream >> h_a[i * ncomps + j];
    }
  }
  return h_a.write();
}

}  // end namespace Omega_h
