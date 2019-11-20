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
#include "Omega_h_for.hpp"
#include "Omega_h_inertia.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

bool is_little_endian_cpu() {
  static std::uint16_t const endian_canary = 0x1;
  std::uint8_t const* p = reinterpret_cast<std::uint8_t const*>(&endian_canary);
  return *p == 0x1;
}

namespace binary {

namespace {

static_assert(sizeof(Int) == 4, "osh format assumes 32 bit Int");
static_assert(sizeof(LO) == 4, "osh format assumes 32 bit LO");
static_assert(sizeof(GO) == 8, "osh format assumes 64 bit GO");
static_assert(sizeof(Real) == 8, "osh format assumes 64 bit Real");

OMEGA_H_INLINE std::uint32_t bswap32(std::uint32_t a) {
#if defined(__GNUC__) && !defined(__CUDA_ARCH__)
  a = __builtin_bswap32(a);
#elif defined(_MSC_VER) && !defined(__CUDA_ARCH__)
  a = _byteswap_ulong(a);
#else
  a = ((a & 0x000000FF) << 24) | ((a & 0x0000FF00) << 8) |
      ((a & 0x00FF0000) >> 8) | ((a & 0xFF000000) >> 24);
#endif
  return a;
}

OMEGA_H_INLINE std::uint64_t bswap64(std::uint64_t a) {
#if defined(__GNUC__) && !defined(__CUDA_ARCH__)
  a = __builtin_bswap64(a);
#elif defined(_MSC_VER) && !defined(__CUDA_ARCH__)
  a = _byteswap_uint64(a);
#else
  a = ((a & 0x00000000000000FFULL) << 56) |
      ((a & 0x000000000000FF00ULL) << 40) |
      ((a & 0x0000000000FF0000ULL) << 24) | ((a & 0x00000000FF000000ULL) << 8) |
      ((a & 0x000000FF00000000ULL) >> 8) | ((a & 0x0000FF0000000000ULL) >> 24) |
      ((a & 0x00FF000000000000ULL) >> 40) | ((a & 0xFF00000000000000ULL) >> 56);
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

unsigned char const magic[2] = {0xa1, 0x1a};

}  // end anonymous namespace

template <typename T>
void swap_bytes(T& ptr) {
  SwapBytes<T>::swap(&ptr);
}

template <typename T>
Read<T> swap_bytes(Read<T> array, bool needs_swapping) {
  if (!needs_swapping) return array;
  Write<T> out = deep_copy(array);
  auto f = OMEGA_H_LAMBDA(LO i) { SwapBytes<T>::swap(&out[i]); };
  parallel_for(out.size(), f, "swap_if_needed");
  return out;
}

template <typename T>
void write_value(std::ostream& stream, T val, bool needs_swapping) {
  if (needs_swapping) swap_bytes(val);
  stream.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
void read_value(std::istream& stream, T& val, bool needs_swapping) {
  stream.read(reinterpret_cast<char*>(&val), sizeof(T));
  if (needs_swapping) swap_bytes(val);
}

template <typename T>
void write_array(std::ostream& stream, Read<T> array, bool is_compressed,
    bool needs_swapping) {
  LO size = array.size();
  write_value(stream, size, needs_swapping);
  Read<T> swapped = swap_bytes(array, needs_swapping);
  HostRead<T> uncompressed(swapped);
  I64 uncompressed_bytes =
      static_cast<I64>(static_cast<std::size_t>(size) * sizeof(T));
#ifdef OMEGA_H_USE_ZLIB
  if (is_compressed) {
    uLong source_bytes = static_cast<uLong>(uncompressed_bytes);
    uLong dest_bytes = ::compressBound(source_bytes);
    auto compressed = new ::Bytef[dest_bytes];
    int ret = ::compress2(compressed, &dest_bytes,
        reinterpret_cast<const ::Bytef*>(nonnull(uncompressed.data())),
        source_bytes, Z_BEST_SPEED);
    OMEGA_H_CHECK(ret == Z_OK);
    I64 compressed_bytes = static_cast<I64>(dest_bytes);
    write_value(stream, compressed_bytes, needs_swapping);
    stream.write(reinterpret_cast<const char*>(compressed), compressed_bytes);
    delete[] compressed;
  } else
#else
  OMEGA_H_CHECK(is_compressed == false);
#endif
  {
    stream.write(reinterpret_cast<const char*>(nonnull(uncompressed.data())),
        uncompressed_bytes);
  }
}

template <typename T>
void read_array(std::istream& stream, Read<T>& array, bool is_compressed,
    bool needs_swapping) {
  LO size;
  read_value(stream, size, needs_swapping);
  OMEGA_H_CHECK(size >= 0);
  I64 uncompressed_bytes =
      static_cast<I64>(static_cast<std::size_t>(size) * sizeof(T));
  HostWrite<T> uncompressed(size);
#ifdef OMEGA_H_USE_ZLIB
  if (is_compressed) {
    I64 compressed_bytes;
    read_value(stream, compressed_bytes, needs_swapping);
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
  array = swap_bytes(Read<T>(uncompressed.write()), needs_swapping);
}

void write(std::ostream& stream, std::string const& val, bool needs_swapping) {
  I32 len = static_cast<I32>(val.length());
  write_value(stream, len, needs_swapping);
  stream.write(val.c_str(), len);
}

void read(std::istream& stream, std::string& val, bool needs_swapping) {
  I32 len;
  read_value(stream, len, needs_swapping);
  OMEGA_H_CHECK(len >= 0);
  val.resize(static_cast<std::size_t>(len));
  stream.read(&val[0], len);
}

static void write_meta(
    std::ostream& stream, Mesh const* mesh, bool needs_swapping) {
  auto family = I8(mesh->family());
  write_value(stream, family, needs_swapping);
  auto dim = I8(mesh->dim());
  write_value(stream, dim, needs_swapping);
  I32 comm_size = mesh->comm()->size();
  write_value(stream, comm_size, needs_swapping);
  I32 comm_rank = mesh->comm()->rank();
  write_value(stream, comm_rank, needs_swapping);
  I8 parting = mesh->parting();
  write_value(stream, parting, needs_swapping);
  I32 nghost_layers = mesh->nghost_layers();
  write_value(stream, nghost_layers, needs_swapping);
  auto hints = mesh->rib_hints();
  I8 have_hints = (hints != nullptr);
  write_value(stream, have_hints, needs_swapping);
  if (have_hints) {
    auto naxes = I32(hints->axes.size());
    write_value(stream, naxes, needs_swapping);
    for (auto axis : hints->axes) {
      for (Int i = 0; i < 3; ++i) write_value(stream, axis[i], needs_swapping);
    }
  }
}

static void read_meta(
    std::istream& stream, Mesh* mesh, Int version, bool needs_swapping) {
  if (version >= 7) {
    I8 family;
    read_value(stream, family, needs_swapping);
    mesh->set_family(Omega_h_Family(family));
  }
  I8 dim;
  read_value(stream, dim, needs_swapping);
  mesh->set_dim(Int(dim));
  I32 comm_size;
  read_value(stream, comm_size, needs_swapping);
  OMEGA_H_CHECK(mesh->comm()->size() == comm_size);
  I32 comm_rank;
  read_value(stream, comm_rank, needs_swapping);
  OMEGA_H_CHECK(mesh->comm()->rank() == comm_rank);
  I8 parting_i8;
  read_value(stream, parting_i8, needs_swapping);
  OMEGA_H_CHECK(parting_i8 == I8(OMEGA_H_ELEM_BASED) ||
                parting_i8 == I8(OMEGA_H_GHOSTED) ||
                parting_i8 == I8(OMEGA_H_VERT_BASED));
  if (version >= 3) {
    I32 nghost_layers;
    read_value(stream, nghost_layers, needs_swapping);
    mesh->set_parting(Omega_h_Parting(parting_i8), nghost_layers, false);
  } else {
    mesh->set_parting(Omega_h_Parting(parting_i8));
  }
  I8 have_hints;
  read_value(stream, have_hints, needs_swapping);
  if (have_hints) {
    I32 naxes;
    read_value(stream, naxes, needs_swapping);
    auto hints = std::make_shared<inertia::Rib>();
    for (I32 i = 0; i < naxes; ++i) {
      Vector<3> axis;
      for (Int j = 0; j < 3; ++j) read_value(stream, axis[j], needs_swapping);
      hints->axes.push_back(axis);
    }
    mesh->set_rib_hints(hints);
  }
  if (version < 6) {
    I8 keeps_canon;
    read_value(stream, keeps_canon, needs_swapping);
  }
}

static void write_tag(std::ostream& stream, TagBase const* tag,
    bool is_compressed, bool needs_swapping) {
  std::string name = tag->name();
  write(stream, name, needs_swapping);
  auto ncomps = I8(tag->ncomps());
  write_value(stream, ncomps, needs_swapping);
  I8 type = tag->type();
  write_value(stream, type, needs_swapping);
  if (is<I8>(tag)) {
    write_array(stream, as<I8>(tag)->array(), is_compressed, needs_swapping);
  } else if (is<I32>(tag)) {
    write_array(stream, as<I32>(tag)->array(), is_compressed, needs_swapping);
  } else if (is<I64>(tag)) {
    write_array(stream, as<I64>(tag)->array(), is_compressed, needs_swapping);
  } else if (is<Real>(tag)) {
    write_array(stream, as<Real>(tag)->array(), is_compressed, needs_swapping);
  } else {
    Omega_h_fail("unexpected tag type in binary write\n");
  }
}

static void read_tag(std::istream& stream, Mesh* mesh, Int d,
    bool is_compressed, I32 version, bool needs_swapping) {
  std::string name;
  read(stream, name, needs_swapping);
  I8 ncomps;
  read_value(stream, ncomps, needs_swapping);
  I8 type;
  read_value(stream, type, needs_swapping);
  if (version < 5) {
    I8 xfer_i8;
    read_value(stream, xfer_i8, needs_swapping);
    if (2 <= version) {
      I8 outflags_i8;
      read_value(stream, outflags_i8, needs_swapping);
    }
  }
  if (type == OMEGA_H_I8) {
    Read<I8> array;
    read_array(stream, array, is_compressed, needs_swapping);
    mesh->add_tag(d, name, ncomps, array, true);
  } else if (type == OMEGA_H_I32) {
    Read<I32> array;
    read_array(stream, array, is_compressed, needs_swapping);
    mesh->add_tag(d, name, ncomps, array, true);
  } else if (type == OMEGA_H_I64) {
    Read<I64> array;
    read_array(stream, array, is_compressed, needs_swapping);
    mesh->add_tag(d, name, ncomps, array, true);
  } else if (type == OMEGA_H_F64) {
    Read<Real> array;
    read_array(stream, array, is_compressed, needs_swapping);
    mesh->add_tag(d, name, ncomps, array, true);
  } else {
    Omega_h_fail("unexpected tag type in binary read\n");
  }
}

static void write_sets(std::ostream& stream, Mesh* mesh, bool needs_swapping) {
  auto n = I32(mesh->class_sets.size());
  write_value(stream, n, needs_swapping);
  for (auto& set : mesh->class_sets) {
    auto& name = set.first;
    write(stream, name, needs_swapping);
    auto npairs = I32(set.second.size());
    write_value(stream, npairs, needs_swapping);
    for (auto& pair : set.second) {
      write_value(stream, pair.dim, needs_swapping);
      write_value(stream, pair.id, needs_swapping);
    }
  }
}

static void read_sets(std::istream& stream, Mesh* mesh, bool needs_swapping) {
  I32 n;
  read_value(stream, n, needs_swapping);
  for (I32 i = 0; i < n; ++i) {
    std::string name;
    read(stream, name, needs_swapping);
    I32 npairs;
    read_value(stream, npairs, needs_swapping);
    for (I32 j = 0; j < npairs; ++j) {
      ClassPair pair;
      read_value(stream, pair.dim, needs_swapping);
      read_value(stream, pair.id, needs_swapping);
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
  bool needs_swapping = !is_little_endian_cpu();
  write_value(stream, is_compressed, needs_swapping);
  write_meta(stream, mesh, needs_swapping);
  LO nverts = mesh->nverts();
  write_value(stream, nverts, needs_swapping);
  for (Int d = 1; d <= mesh->dim(); ++d) {
    auto down = mesh->ask_down(d, d - 1);
    write_array(stream, down.ab2b, is_compressed, needs_swapping);
    if (d > 1) {
      write_array(stream, down.codes, is_compressed, needs_swapping);
    }
  }
  for (Int d = 0; d <= mesh->dim(); ++d) {
    auto nsaved_tags = mesh->ntags(d);
    write_value(stream, nsaved_tags, needs_swapping);
    for (Int i = 0; i < mesh->ntags(d); ++i) {
      write_tag(stream, mesh->get_tag(d, i), is_compressed, needs_swapping);
    }
    if (mesh->comm()->size() > 1) {
      auto owners = mesh->ask_owners(d);
      write_array(stream, owners.ranks, is_compressed, needs_swapping);
      write_array(stream, owners.idxs, is_compressed, needs_swapping);
    }
  }
  write_sets(stream, mesh, needs_swapping);
  I8 has_parents = mesh->has_any_parents();
  write_value(stream, has_parents, needs_swapping);
  if (has_parents) {
    for (Int d = 0; d <= mesh->dim(); ++d) {
      auto parents = mesh->ask_parents(d);
      write_array(stream, parents.parent_idx, is_compressed, needs_swapping);
      write_array(stream, parents.codes, is_compressed, needs_swapping);
    }
  }
  end_code();
}

void read(std::istream& stream, Mesh* mesh, I32 version) {
  ScopedTimer timer("binary::read(istream, mesh, version)");
  unsigned char magic_in[2];
  stream.read(reinterpret_cast<char*>(magic_in), sizeof(magic));
  OMEGA_H_CHECK(magic_in[0] == magic[0]);
  OMEGA_H_CHECK(magic_in[1] == magic[1]);
  bool needs_swapping = !is_little_endian_cpu();
  if (version == -1) read_value(stream, version, needs_swapping);
  OMEGA_H_CHECK(version >= 1);
  OMEGA_H_CHECK(version <= latest_version);
  I8 is_compressed;
  read_value(stream, is_compressed, needs_swapping);
#ifndef OMEGA_H_USE_ZLIB
  OMEGA_H_CHECK(!is_compressed);
#endif
  read_meta(stream, mesh, version, needs_swapping);
  LO nverts;
  read_value(stream, nverts, needs_swapping);
  mesh->set_verts(nverts);
  for (Int d = 1; d <= mesh->dim(); ++d) {
    Adj down;
    read_array(stream, down.ab2b, is_compressed, needs_swapping);
    if (d > 1) {
      read_array(stream, down.codes, is_compressed, needs_swapping);
    }
    mesh->set_ents(d, down);
  }
  for (Int d = 0; d <= mesh->dim(); ++d) {
    Int ntags;
    read_value(stream, ntags, needs_swapping);
    for (Int i = 0; i < ntags; ++i) {
      read_tag(stream, mesh, d, is_compressed, version, needs_swapping);
    }
    if (mesh->comm()->size() > 1) {
      Remotes owners;
      read_array(stream, owners.ranks, is_compressed, needs_swapping);
      read_array(stream, owners.idxs, is_compressed, needs_swapping);
      mesh->set_owners(d, owners);
    }
  }
  if (version >= 8) {
    read_sets(stream, mesh, needs_swapping);
  }
  if (version >= 9) {
    I8 has_parents;
    read_value(stream, has_parents, needs_swapping);
    if (has_parents) {
      for (Int d = 0; d <= mesh->dim(); ++d) {
        Parents parents;
        read_array(stream, parents.parent_idx, is_compressed, needs_swapping);
        read_array(stream, parents.codes, is_compressed, needs_swapping);
        mesh->set_parents(d, parents);
      }
    }
  }
}

static void write_int_file(
    filesystem::path const& filepath, Mesh* mesh, I32 value) {
  if (mesh->comm()->rank() == 0) {
    std::ofstream file(filepath.c_str());
    OMEGA_H_CHECK(file.is_open());
    file << value << '\n';
  }
}

static void write_nparts(filesystem::path const& path, Mesh* mesh) {
  write_int_file(path / "nparts", mesh, mesh->comm()->size());
}

static void write_version(filesystem::path const& path, Mesh* mesh) {
  write_int_file(path / "version", mesh, latest_version);
}

I32 read_nparts(filesystem::path const& path, CommPtr comm) {
  I32 nparts;
  if (comm->rank() == 0) {
    auto const filepath = path / "nparts";
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

I32 read_version(filesystem::path const& path, CommPtr comm) {
  I32 version;
  if (comm->rank() == 0) {
    auto const filepath = path / "version";
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

void write(filesystem::path const& path, Mesh* mesh) {
  begin_code("binary::write(path,Mesh)");
  if (path.extension().string() != ".osh" && can_print(mesh)) {
    std::cout
        << "it is strongly recommended to end Omega_h paths in \".osh\",\n";
    std::cout << "instead of just \"" << path << "\"\n";
  }
  filesystem::create_directory(path);
  mesh->comm()->barrier();
  auto filepath = path;
  filepath /= std::to_string(mesh->comm()->rank());
  filepath += ".osh";
  std::ofstream file(filepath.c_str());
  OMEGA_H_CHECK(file.is_open());
  write(file, mesh);
  write_nparts(path, mesh);
  write_version(path, mesh);
  mesh->comm()->barrier();
  end_code();
}

void read_in_comm(
    filesystem::path const& path, CommPtr comm, Mesh* mesh, I32 version) {
  ScopedTimer timer("binary::read_in_comm(path, comm, mesh, version)");
  mesh->set_comm(comm);
  auto filepath = path;
  filepath /= std::to_string(mesh->comm()->rank());
  if (version != -1) filepath += ".osh";
  std::ifstream file(filepath.c_str(), std::ios::binary);
  OMEGA_H_CHECK(file.is_open());
  read(file, mesh, version);
}

I32 read(filesystem::path const& path, CommPtr comm, Mesh* mesh, bool strict) {
  ScopedTimer timer("binary::read(path, comm, mesh, strict)");
  auto const nparts = read_nparts(path, comm);
  auto const version = read_version(path, comm);
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
    auto const in_subcomm = (comm->rank() < nparts);
    auto const subcomm = comm->split(I32(!in_subcomm), 0);
    if (in_subcomm) {
      read_in_comm(path, subcomm, mesh, version);
    }
    mesh->set_comm(comm);
  }
  return nparts;
}

Mesh read(filesystem::path const& path, Library* lib, bool strict) {
  ScopedTimer timer("binary::read(path, lib, strict)");
  return binary::read(path, lib->world(), strict);
}

Mesh read(filesystem::path const& path, CommPtr comm, bool strict) {
  ScopedTimer timer("binary::read(path, comm, strict)");
  auto mesh = Mesh(comm->library());
  binary::read(path, comm, &mesh, strict);
  return mesh;
}

#define OMEGA_H_INST(T)                                                        \
  template void swap_bytes(T&);                                                \
  template Read<T> swap_bytes(Read<T> array, bool is_little_endian);           \
  template void write_value(std::ostream& stream, T val, bool);                \
  template void read_value(std::istream& stream, T& val, bool);                \
  template void write_array(std::ostream& stream, Read<T> array, bool, bool);  \
  template void read_array(                                                    \
      std::istream& stream, Read<T>& array, bool is_compressed, bool);
OMEGA_H_INST(I8)
OMEGA_H_INST(I32)
OMEGA_H_INST(I64)
OMEGA_H_INST(Real)
#undef OMEGA_H_INST

// for VTK compression headers
template void swap_bytes(std::uint64_t&);

}  // end namespace binary

void write_reals_txt(filesystem::path const& filename, Reals a, Int ncomps) {
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

Reals read_reals_txt(filesystem::path const& filename, LO n, Int ncomps) {
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

OMEGA_H_DLL Mesh read_mesh_file(filesystem::path const& path, CommPtr comm) {
  auto const extension = path.extension().string();
  if (extension == ".osh") {
    return binary::read(path, comm);
  } else if (extension == ".meshb") {
#ifdef OMEGA_H_USE_LIBMESHB
    Mesh mesh(comm->library());
    meshb::read(&mesh, path);
    mesh.set_comm(comm);
    return mesh;
#else
    Omega_h_fail(
        "Omega_h: Can't read %s without reconfiguring with "
        "OMEGA_H_USE_libMeshb=ON\n",
        path.c_str());
    OMEGA_H_NORETURN(Mesh());
#endif
  } else if (extension == ".exo" || extension == ".e" || extension == ".g") {
#ifdef OMEGA_H_USE_SEACASEXODUS
    Mesh mesh(comm->library());
    auto file = exodus::open(path);
    exodus::read_mesh(file, &mesh);
    mesh.set_comm(comm);
    return mesh;
#else
    Omega_h_fail(
        "Omega_h: Can't read %s without reconfiguring with "
        "OMEGA_H_USE_SEACASExodus=ON\n",
        path.c_str());
    OMEGA_H_NORETURN(Mesh());
#endif
  } else if (extension == ".msh") {
    return gmsh::read(path, comm);
  } else if (extension == ".pvtu") {
    Mesh mesh(comm->library());
    vtk::read_parallel(path, comm, &mesh);
    mesh.set_comm(comm);
    return mesh;
  } else if (extension == ".vtu") {
    Mesh mesh(comm->library());
    std::ifstream stream(path.c_str());
    OMEGA_H_CHECK(stream.is_open());
    vtk::read_vtu(stream, comm, &mesh);
    return mesh;
  } else {
    Omega_h_fail("Unknown file extension \"%s\" on \"%s\"\n", extension.c_str(),
        path.c_str());
  }
}

}  // end namespace Omega_h
