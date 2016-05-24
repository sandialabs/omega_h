bool is_little_endian_cpu()
{
  static std::uint16_t const endian_canary = 0x1;
  std::uint8_t const* p = reinterpret_cast<std::uint8_t const*>(&endian_canary);
  return *p == 0x1;
}

namespace file {

namespace {

static_assert(sizeof(LO) == 4, "osh format assumes 32 bit LO");
static_assert(sizeof(GO) == 8, "osh format assumes 64 bit GO");
static_assert(sizeof(Real) == 8, "osh format assumes 64 bit Real");

enum {
  OSH_I8  = 0,
  OSH_I32 = 2,
  OSH_I64 = 3,
  OSH_F64 = 5,
};

INLINE std::uint32_t bswap32(std::uint32_t a)
{
#ifdef USE_CUDA
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
#ifdef USE_CUDA
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

} //end anonymous namespace

template <typename T>
void write_binary(std::ostream& stream, T val) {
  swap_if_needed(val);
  stream.write(reinterpret_cast<const char*>(&val), sizeof(T));
}

template <typename T>
void read_binary(std::istream& stream, T& val) {
  stream.read(reinterpret_cast<char*>(&val), sizeof(T));
  swap_if_needed(val);
}

template <typename T>
void write_binary(std::ostream& stream, Read<T> array) {
  LO size = array.size();
  write_binary(stream, size);
  Read<T> swapped = swap_if_needed(array);
  HostRead<T> uncompressed(swapped);
  I64 uncompressed_bytes = static_cast<I64>(
      static_cast<std::size_t>(size) * sizeof(T));
#ifdef USE_ZLIB
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
  write_binary(stream, compressed_bytes);
  stream.write(reinterpret_cast<const char*>(compressed), compressed_bytes);
  delete [] compressed;
#else
  stream.write(reinterpret_cast<const char*>(&uncompressed[0]),
      uncompressed_bytes);
#endif
}

template <typename T>
void read_binary(std::istream& stream, Read<T>& array,
    bool is_compressed) {
  LO size;
  read_binary(stream, size);
  I64 uncompressed_bytes = static_cast<I64>(
      static_cast<std::size_t>(size) * sizeof(T));
  HostWrite<T> uncompressed(size);
#ifdef USE_ZLIB
  if (is_compressed) {
    I64 compressed_bytes;
    read_binary(stream, compressed_bytes);
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
#endif
  {
    stream.read(reinterpret_cast<char*>(&uncompressed[0]),
        uncompressed_bytes);
  }
  array = swap_if_needed(Read<T>(uncompressed.write()));
}

#define INST_T(T) \
template void write_binary(std::ostream& stream, T val); \
template void read_binary(std::istream& stream, T& val); \
template void write_binary(std::ostream& stream, Read<T> array); \
template void read_binary(std::istream& stream, Read<T>& array, \
    bool is_compressed);
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T

} //end namespace file
