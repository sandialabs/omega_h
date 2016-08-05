#include "file.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <fstream>
#include <iostream>

#ifdef OSH_USE_ZLIB
#include <zlib.h>
#endif

#include "array.hpp"
#include "inertia.hpp"
#include "loop.hpp"
#include "tag.hpp"

namespace osh {

bool is_little_endian_cpu() {
  static std::uint16_t const endian_canary = 0x1;
  std::uint8_t const* p = reinterpret_cast<std::uint8_t const*>(&endian_canary);
  return *p == 0x1;
}

void safe_mkdir(const char* path) {
  mode_t const mode = S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST) {
    osh_fail("omega_h could not create directory \"%s\"\n", path);
  }
}

bool directory_exists(const char* path) {
  struct stat info;
  if (stat(path, &info) != 0) return false;
  CHECK(info.st_mode & S_IFDIR);
  return true;
}

std::string parent_path(std::string const& path) {
  auto pos = path.find_last_of('/');
  CHECK(pos != std::string::npos);
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

INLINE std::uint32_t bswap32(std::uint32_t a) {
#ifdef OSH_USE_CUDA
  a = ((a & 0x000000FF) << 24) | ((a & 0x0000FF00) << 8) |
      ((a & 0x00FF0000) >> 8) | ((a & 0xFF000000) >> 24);
#else
  a = __builtin_bswap32(a);
#endif
  return a;
}

INLINE std::uint64_t bswap64(std::uint64_t a) {
#ifdef OSH_USE_CUDA
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
  INLINE static void swap(T*) {}
};

template <typename T>
struct SwapBytes<T, 4> {
  INLINE static void swap(T* ptr) {
    std::uint32_t* p2 = reinterpret_cast<std::uint32_t*>(ptr);
    *p2 = bswap32(*p2);
  }
};

template <typename T>
struct SwapBytes<T, 8> {
  INLINE static void swap(T* ptr) {
    std::uint64_t* p2 = reinterpret_cast<std::uint64_t*>(ptr);
    *p2 = bswap64(*p2);
  }
};

template <typename T>
INLINE void swap_bytes(T* ptr) {
  SwapBytes<T>::swap(ptr);
}

unsigned char const magic[2] = {0xa1, 0x1a};
I32 latest_version = 1;

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
  auto f = LAMBDA(LO i) { swap_bytes(&out[i]); };
  parallel_for(out.size(), f);
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
#ifdef OSH_USE_ZLIB
  uLong source_bytes = static_cast<uLong>(uncompressed_bytes);
  uLong dest_bytes = ::compressBound(source_bytes);
  auto compressed = new Bytef[dest_bytes];
  int ret = ::compress2(compressed, &dest_bytes,
                        reinterpret_cast<const Bytef*>(uncompressed.data()),
                        source_bytes, Z_BEST_SPEED);
  CHECK(ret == Z_OK);
  I64 compressed_bytes = static_cast<I64>(dest_bytes);
  write_value(stream, compressed_bytes);
  stream.write(reinterpret_cast<const char*>(compressed), compressed_bytes);
  delete[] compressed;
#else
  stream.write(reinterpret_cast<const char*>(uncompressed.data()),
               uncompressed_bytes);
#endif
}

template <typename T>
void read_array(std::istream& stream, Read<T>& array, bool is_compressed) {
  LO size;
  read_value(stream, size);
  CHECK(size >= 0);
  I64 uncompressed_bytes =
      static_cast<I64>(static_cast<std::size_t>(size) * sizeof(T));
  HostWrite<T> uncompressed(size);
#ifdef OSH_USE_ZLIB
  if (is_compressed) {
    I64 compressed_bytes;
    read_value(stream, compressed_bytes);
    CHECK(compressed_bytes >= 0);
    auto compressed = new Bytef[compressed_bytes];
    stream.read(reinterpret_cast<char*>(compressed), compressed_bytes);
    uLong dest_bytes = static_cast<uLong>(uncompressed_bytes);
    uLong source_bytes = static_cast<uLong>(compressed_bytes);
    Bytef* uncompressed_ptr = reinterpret_cast<Bytef*>(uncompressed.data());
    int ret =
        ::uncompress(uncompressed_ptr, &dest_bytes, compressed, source_bytes);
    CHECK(ret == Z_OK);
    CHECK(dest_bytes == static_cast<uLong>(uncompressed_bytes));
    delete[] compressed;
  } else
#else
  CHECK(is_compressed == false);
#endif
  {
    stream.read(reinterpret_cast<char*>(uncompressed.data()),
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
  CHECK(len >= 0);
  val.resize(static_cast<std::size_t>(len));
  stream.read(&val[0], len);
}

static void write_meta(std::ostream& stream, Mesh const* mesh) {
  auto dim = I8(mesh->dim());
  write_value(stream, dim);
  I32 comm_size = mesh->comm()->size();
  write_value(stream, comm_size);
  I32 comm_rank = mesh->comm()->rank();
  write_value(stream, comm_rank);
  I8 parting = mesh->parting();
  write_value(stream, parting);
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
  I8 keeps_canon = mesh->keeps_canonical_globals();
  write_value(stream, keeps_canon);
}

static void read_meta(std::istream& stream, Mesh* mesh) {
  I8 dim;
  read_value(stream, dim);
  mesh->set_dim(Int(dim));
  I32 comm_size;
  read_value(stream, comm_size);
  CHECK(mesh->comm()->size() == comm_size);
  I32 comm_rank;
  read_value(stream, comm_rank);
  CHECK(mesh->comm()->rank() == comm_rank);
  I8 parting;
  read_value(stream, parting);
  mesh->set_parting(osh_parting(parting));
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
  I8 keeps_canon;
  read_value(stream, keeps_canon);
  mesh->keep_canonical_globals(bool(keeps_canon));
}

static void write_tag(std::ostream& stream, TagBase const* tag) {
  std::string name = tag->name();
  write(stream, name);
  auto ncomps = I8(tag->ncomps());
  write_value(stream, ncomps);
  I8 type = tag->type();
  write_value(stream, type);
  I8 xfer_int = static_cast<I8>(tag->xfer());
  write_value(stream, xfer_int);
  if (is<I8>(tag)) {
    write_array(stream, to<I8>(tag)->array());
  } else if (is<I32>(tag)) {
    write_array(stream, to<I32>(tag)->array());
  } else if (is<I64>(tag)) {
    write_array(stream, to<I64>(tag)->array());
  } else if (is<Real>(tag)) {
    write_array(stream, to<Real>(tag)->array());
  } else {
    osh_fail("unexpected tag type in binary write\n");
  }
}

static void read_tag(std::istream& stream, Mesh* mesh, Int d,
                     bool is_compressed) {
  std::string name;
  read(stream, name);
  I8 ncomps;
  read_value(stream, ncomps);
  I8 type;
  read_value(stream, type);
  I8 xfer_i8;
  read_value(stream, xfer_i8);
  Int xfer = static_cast<Int>(xfer_i8);
  if (type == OSH_I8) {
    Read<I8> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, xfer, array);
  } else if (type == OSH_I32) {
    Read<I32> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, xfer, array);
  } else if (type == OSH_I64) {
    Read<I64> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, xfer, array);
  } else if (type == OSH_F64) {
    Read<Real> array;
    read_array(stream, array, is_compressed);
    mesh->add_tag(d, name, ncomps, xfer, array);
  } else {
    osh_fail("unexpected tag type in binary read\n");
  }
}

void write(std::ostream& stream, Mesh* mesh) {
  stream.write(reinterpret_cast<const char*>(magic), sizeof(magic));
  write_value(stream, latest_version);
#ifdef OSH_USE_ZLIB
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
    Int ntags = mesh->ntags(d);
    write_value(stream, ntags);
    for (Int i = 0; i < ntags; ++i) {
      write_tag(stream, mesh->get_tag(d, i));
    }
    if (mesh->comm()->size() > 1) {
      auto owners = mesh->ask_owners(d);
      write_array(stream, owners.ranks);
      write_array(stream, owners.idxs);
    }
  }
}

void read(std::istream& stream, Mesh* mesh) {
  unsigned char magic_in[2];
  stream.read(reinterpret_cast<char*>(magic_in), sizeof(magic));
  CHECK(magic_in[0] == magic[0]);
  CHECK(magic_in[1] == magic[1]);
  I32 version;
  read_value(stream, version);
  CHECK(version >= 1);
  CHECK(version <= latest_version);
  I8 is_compressed;
  read_value(stream, is_compressed);
#ifndef OSH_USE_ZLIB
  CHECK(!is_compressed);
#endif
  read_meta(stream, mesh);
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
      read_tag(stream, mesh, d, is_compressed);
    }
    if (mesh->comm()->size() > 1) {
      Remotes owners;
      read_array(stream, owners.ranks, is_compressed);
      read_array(stream, owners.idxs, is_compressed);
      mesh->set_owners(d, owners);
    }
  }
}

static void write_nparts(std::string const& path, Mesh* mesh) {
  auto filepath = path + "/nparts";
  std::ofstream file(filepath.c_str());
  CHECK(file.is_open());
  file << mesh->comm()->size() << '\n';
}

static I32 read_nparts(std::string const& path) {
  auto filepath = path + "/nparts";
  std::ifstream file(filepath.c_str());
  if (!file.is_open()) {
    osh_fail("could not open file \"%s\"\n", filepath.c_str());
  }
  I32 nparts;
  file >> nparts;
  return nparts;
}

static bool ends_with(std::string const& path, std::string const& ext) {
  if (path.length() < ext.length()) return false;
  return path.substr(path.length() - ext.length(), ext.length()) == ext;
}

void write(std::string const& path, Mesh* mesh) {
  if (!ends_with(path, ".osh") && mesh->comm()->rank() == 0) {
    std::cout
        << "it is strongly recommended to end Omega_h paths in \".osh\",\n";
    std::cout << "instead of just \"" << path << "\"\n";
  }
  safe_mkdir(path.c_str());
  mesh->comm()->barrier();
  auto filepath = path + "/" + to_string(mesh->comm()->rank());
  std::ofstream file(filepath.c_str());
  CHECK(file.is_open());
  write(file, mesh);
  write_nparts(path, mesh);
}

static void read2(std::string const& path, CommPtr comm, Mesh* mesh) {
  mesh->set_comm(comm);
  auto filepath = path + "/" + to_string(mesh->comm()->rank());
  std::ifstream file(filepath.c_str());
  CHECK(file.is_open());
  read(file, mesh);
}

void read(std::string const& path, CommPtr comm, Mesh* mesh) {
  auto nparts = read_nparts(path);
  if (nparts > comm->size()) {
    osh_fail(
        "path \"%s\" contains %d parts, but only %d ranks are reading it\n",
        path.c_str(), nparts, comm->size());
  }
  bool in_subcomm = (comm->rank() < nparts);
  auto subcomm = comm->split(I32(!in_subcomm), 0);
  if (in_subcomm) {
    read2(path, subcomm, mesh);
  }
  mesh->set_comm(comm);
}

#define INST(T)                                                          \
  template void swap_if_needed(T& val, bool is_little_endian);           \
  template Read<T> swap_if_needed(Read<T> array, bool is_little_endian); \
  template void write_value(std::ostream& stream, T val);                \
  template void read_value(std::istream& stream, T& val);                \
  template void write_array(std::ostream& stream, Read<T> array);        \
  template void read_array(std::istream& stream, Read<T>& array,         \
                           bool is_compressed);
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

// for VTK compression headers
template void swap_if_needed(std::size_t& val, bool is_little_endian);

}  // end namespace binary

}  // end namespace osh
