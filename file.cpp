#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */

bool is_little_endian_cpu()
{
  static std::uint16_t const endian_canary = 0x1;
  std::uint8_t const* p = reinterpret_cast<std::uint8_t const*>(&endian_canary);
  return *p == 0x1;
}

void safe_mkdir(const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST) {
    fail("omega_h could not create directory \"%s\"\n", path);
  }
}

std::string parent_path(std::string const& path) {
  auto pos = path.find_last_of('/');
  CHECK(pos != std::string::npos);
  return path.substr(0, pos);
}

std::string path_leaf_name(std::string const& path) {
  auto pos = path.find_last_of('/');
  if (pos == std::string::npos)
    return path;
  return path.substr(pos + 1, std::string::npos);
}

namespace binary {

namespace {

static_assert(sizeof(Int) == 4, "osh format assumes 32 bit Int");
static_assert(sizeof(LO) == 4, "osh format assumes 32 bit LO");
static_assert(sizeof(GO) == 8, "osh format assumes 64 bit GO");
static_assert(sizeof(Real) == 8, "osh format assumes 64 bit Real");

INLINE std::uint32_t bswap32(std::uint32_t a)
{
#ifdef OSH_USE_CUDA
  a = ((a & 0x000000FF) << 24) |
      ((a & 0x0000FF00) <<  8) |
      ((a & 0x00FF0000) >>  8) |
      ((a & 0xFF000000) >> 24);
#else
  a = __builtin_bswap32(a);
#endif
  return a;
}

INLINE std::uint64_t bswap64(std::uint64_t a)
{
#ifdef OSH_USE_CUDA
  a = ((a & 0x00000000000000FFULL) << 56) |
      ((a & 0x000000000000FF00ULL) << 40) |
      ((a & 0x0000000000FF0000ULL) << 24) |
      ((a & 0x00000000FF000000ULL) <<  8) |
      ((a & 0x000000FF00000000ULL) >>  8) |
      ((a & 0x0000FF0000000000ULL) >> 24) |
      ((a & 0x00FF000000000000ULL) >> 40) |
      ((a & 0xFF00000000000000ULL) >> 56);
#else
  a = __builtin_bswap64(a);
#endif
  return a;
}

template <typename T, size_t size = sizeof(T)>
struct SwapBytes;

template <typename T>
struct SwapBytes<T, 1> {
  INLINE static void swap(T*) {
  }
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

template <typename T>
void swap_if_needed(T& val) {
  if (!is_little_endian_cpu()) {
    swap_bytes(&val);
  }
}

template <typename T>
static Read<T> swap_if_needed(Read<T> array) {
  if (is_little_endian_cpu()) {
    return array;
  }
  Write<T> out = deep_copy(array);
  auto f = LAMBDA(LO i) {
    swap_bytes(&out[i]);
  };
  parallel_for(out.size(), f);
  return out;
}

unsigned char const magic[2] = {0xa1,0x1a};
Int latest_version = 1;

} //end anonymous namespace

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
  Read<T> swapped = swap_if_needed(array);
  HostRead<T> uncompressed(swapped);
  I64 uncompressed_bytes = static_cast<I64>(
      static_cast<std::size_t>(size) * sizeof(T));
#ifdef OSH_USE_ZLIB
  uLong source_bytes = static_cast<uLong>(uncompressed_bytes);
  uLong dest_bytes = ::compressBound(source_bytes);
  Bytef* compressed = new Bytef[dest_bytes];
  int ret = ::compress2(
      compressed,
      &dest_bytes,
      reinterpret_cast<const Bytef*>(&uncompressed[0]),
      source_bytes,
      9);
  CHECK(ret == Z_OK);
  I64 compressed_bytes = static_cast<I64>(dest_bytes);
  write_value(stream, compressed_bytes);
  stream.write(reinterpret_cast<const char*>(compressed), compressed_bytes);
  delete [] compressed;
#else
  stream.write(reinterpret_cast<const char*>(&uncompressed[0]),
      uncompressed_bytes);
#endif
}

template <typename T>
void read_array(std::istream& stream, Read<T>& array,
    bool is_compressed) {
  LO size;
  read_value(stream, size);
  CHECK(size >= 0);
  I64 uncompressed_bytes = static_cast<I64>(
      static_cast<std::size_t>(size) * sizeof(T));
  HostWrite<T> uncompressed(size);
#ifdef OSH_USE_ZLIB
  if (is_compressed) {
    I64 compressed_bytes;
    read_value(stream, compressed_bytes);
    CHECK(compressed_bytes >= 0);
    Bytef* compressed = new Bytef[compressed_bytes];
    stream.read(reinterpret_cast<char*>(compressed), compressed_bytes);
    uLong dest_bytes = static_cast<uLong>(uncompressed_bytes);
    uLong source_bytes = static_cast<uLong>(compressed_bytes);
    int ret = ::uncompress(
        reinterpret_cast<Bytef*>(&uncompressed[0]),
        &dest_bytes,
        compressed,
        source_bytes);
    CHECK(ret == Z_OK);
    CHECK(dest_bytes == static_cast<uLong>(uncompressed_bytes));
    delete [] compressed;
  } else
#else
  CHECK(is_compressed == false);
#endif
  {
    stream.read(reinterpret_cast<char*>(&uncompressed[0]),
        uncompressed_bytes);
  }
  array = swap_if_needed(Read<T>(uncompressed.write()));
}

void write(std::ostream& stream, std::string const& val)
{
  I32 len = static_cast<I32>(val.length());
  write_value(stream, len);
  stream.write(val.c_str(), len);
}

void read(std::istream& stream, std::string& val)
{
  I32 len;
  read_value(stream, len);
  CHECK(len >= 0);
  val.resize(static_cast<std::size_t>(len));
  stream.read(&val[0], len);
}

void write(std::ostream& stream, Mesh& mesh) {
  stream.write(reinterpret_cast<const char*>(magic), sizeof(magic));
  write_value(stream, latest_version);
#ifdef OSH_USE_ZLIB
  I8 is_compressed = true;
#else
  I8 is_compressed = false;
#endif
  write_value(stream, is_compressed);
  Int dim = mesh.dim();
  write_value(stream, dim);
  I32 comm_size = mesh.comm()->size();
  write_value(stream, comm_size);
  Int partition = mesh.partition();
  write_value(stream, partition);
  LO nverts = mesh.nverts();
  write_value(stream, nverts);
  for (Int d = 1; d <= dim; ++d) {
    LOs down = mesh.ask_down(d, d - 1).ab2b;
    write_array(stream, down);
  }
  for (Int d = 0; d <= dim; ++d) {
    Int ntags = mesh.ntags(d);
    write_value(stream, ntags);
    for (Int i = 0; i < ntags; ++i) {
      auto tag = mesh.get_tag(d, i);
      write(stream, tag->name());
      Int ncomps = tag->ncomps();
      write_value(stream, ncomps);
      write_value(stream, static_cast<Int>(tag->type()));
      write_value(stream, static_cast<Int>(tag->xfer()));
      if (is<I8>(tag)) {
        write_array(stream, to<I8>(tag)->array());
      } else if (is<I32>(tag)) {
        write_array(stream, to<I32>(tag)->array());
      } else if (is<I64>(tag)) {
        write_array(stream, to<I64>(tag)->array());
      } else if (is<Real>(tag)) {
        write_array(stream, to<Real>(tag)->array());
      } else {
        fail("unexpected tag type in binary write\n");
      }
    }
    if (mesh.comm()->size() > 1) {
      auto owners = mesh.ask_owners(dim);
      write_array(stream, owners.ranks);
      write_array(stream, owners.idxs);
    }
  }
}

void read(std::istream& stream, Mesh& mesh) {
  unsigned char magic_in[2];
  stream.read(reinterpret_cast<char*>(magic_in), sizeof(magic));
  CHECK(magic_in[0] == magic[0]);
  CHECK(magic_in[1] == magic[1]);
  Int version;
  read_value(stream, version);
  I8 is_compressed;
  read_value(stream, is_compressed);
#ifdef OSH_USE_ZLIB
  CHECK(!is_compressed);
#endif
  CHECK(version >= 1);
  CHECK(version <= latest_version);
  Int dim;
  read_value(stream, dim);
  mesh.set_dim(dim);
  I32 comm_size;
  read_value(stream, comm_size);
  CHECK(comm_size == mesh.comm()->size());
  Int partition;
  read_value(stream, partition);
  mesh.set_partition(static_cast<Partition>(partition));
  LO nverts;
  read_value(stream, nverts);
  for (Int d = 1; d <= dim; ++d) {
    LOs down;
    read_array(stream, down, is_compressed);
    mesh.set_ents(d, down);
  }
  for (Int d = 0; d <= dim; ++d) {
    Int ntags;
    read_value(stream, ntags);
    for (Int i = 0; i < ntags; ++i) {
      std::string name;
      read(stream, name);
      Int ncomps;
      read_value(stream, ncomps);
      Int type;
      read_value(stream, type);
      Int xfer_int;
      read_value(stream, xfer_int);
      Xfer xfer = static_cast<Xfer>(xfer_int);
      if (type == OSH_I8) {
        Read<I8> array;
        read_array(stream, array, is_compressed);
        mesh.add_tag(d, name, ncomps, xfer, array);
      } else if (type == OSH_I32) {
        Read<I32> array;
        read_array(stream, array, is_compressed);
        mesh.add_tag(d, name, ncomps, xfer, array);
      } else if (type == OSH_I64) {
        Read<I64> array;
        read_array(stream, array, is_compressed);
        mesh.add_tag(d, name, ncomps, xfer, array);
      } else if (type == OSH_F64) {
        Read<Real> array;
        read_array(stream, array, is_compressed);
        mesh.add_tag(d, name, ncomps, xfer, array);
      } else {
        fail("unexpected tag type in binary write\n");
      }
    }
    if (mesh.comm()->size() > 1) {
      Remotes owners;
      read_array(stream, owners.ranks, is_compressed);
      read_array(stream, owners.idxs, is_compressed);
      mesh.set_owners(d, owners);
    }
  }
}

#define INST_T(T) \
template void write_value(std::ostream& stream, T val); \
template void read_value(std::istream& stream, T& val); \
template void write_array(std::ostream& stream, Read<T> array); \
template void read_array(std::istream& stream, Read<T>& array, \
    bool is_compressed);
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T

} //end namespace file
