#ifndef OMEGA_H_FILE_HPP
#define OMEGA_H_FILE_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include <Omega_h_config.h>
#include <Omega_h_defines.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_comm.hpp>

namespace Omega_h {

class Mesh;

bool ends_with(std::string const& s, std::string const& suffix);

#ifdef OMEGA_H_USE_LIBMESHB
namespace meshb {
void read(Mesh* mesh, const char* filepath);
void write(Mesh* mesh, const char* filepath, int version);
}  // namespace meshb
#endif

#ifdef OMEGA_H_USE_SEACASEXODUS
namespace exodus {
enum ClassifyWith {
  NODE_SETS = 0x1,
  SIDE_SETS = 0x2,
};
void read(std::string const& path, Mesh* mesh, bool verbose = false,
    int classify_with = NODE_SETS | SIDE_SETS);
void write(std::string const& path, Mesh* mesh, bool verbose = false,
    int classify_with = NODE_SETS | SIDE_SETS);
}  // namespace exodus
#endif

namespace gmsh {
void read(std::istream& stream, Mesh* mesh);
void read(std::string const& filename, Mesh* mesh);
}  // namespace gmsh

namespace vtk {
void write_vtu(std::ostream& stream, Mesh* mesh, Int cell_dim);
void write_vtu(std::string const& filename, Mesh* mesh, Int cell_dim);
void write_parallel(std::string const& path, Mesh* mesh, Int cell_dim);
class Writer {
  Mesh* mesh_;
  std::string root_path_;
  Int cell_dim_;
  Int step_;
  std::streampos pvd_pos_;

 public:
  Writer();
  Writer(Writer const&);
  Writer& operator=(Writer const&);
  ~Writer();
  Writer(Mesh* mesh, std::string const& root_path, Int cell_dim);
  void write(Real time);
  void write();
};
class FullWriter {
  std::vector<Writer> writers_;

 public:
  FullWriter();
  FullWriter(FullWriter const&);
  FullWriter& operator=(FullWriter const&);
  ~FullWriter();
  FullWriter(Mesh* mesh, std::string const& root_path);
  void write(Real time);
  void write();
};
}  // end namespace vtk

bool is_little_endian_cpu();
void safe_mkdir(const char* path);
bool directory_exists(const char* path);
std::string parent_path(std::string const& path);
std::string path_leaf_name(std::string const& path);

namespace binary {

void write(std::string const& path, Mesh* mesh);
I32 read(std::string const& path, CommPtr comm, Mesh* mesh);
I32 read_nparts(std::string const& path, CommPtr comm);
I32 read_version(std::string const& path, CommPtr comm);
void read_in_comm(
    std::string const& path, CommPtr comm, Mesh* mesh, I32 version);

constexpr I32 latest_version = 4;

template <typename T>
void swap_if_needed(T& val, bool is_little_endian = true);
template <typename T>
Read<T> swap_if_needed(Read<T> array, bool is_little_endian);

template <typename T>
void write_value(std::ostream& stream, T val);
template <typename T>
void read_value(std::istream& stream, T& val);
template <typename T>
void write_array(std::ostream& stream, Read<T> array);
template <typename T>
void read_array(std::istream& stream, Read<T>& array, bool is_compressed);

void write(std::ostream& stream, std::string const& val);
void read(std::istream& stream, std::string& val);

void write(std::ostream& stream, Mesh* mesh);
void read(std::istream& stream, Mesh* mesh, I32 version);

#define INST_DECL(T)                                                           \
  extern template void swap_if_needed(T& val, bool is_little_endian);          \
  extern template Read<T> swap_if_needed(                                      \
      Read<T> array, bool is_little_endian);                                   \
  extern template void write_value(std::ostream& stream, T val);               \
  extern template void read_value(std::istream& stream, T& val);               \
  extern template void write_array(std::ostream& stream, Read<T> array);       \
  extern template void read_array(                                             \
      std::istream& stream, Read<T>& array, bool is_compressed);
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL
// for VTK compression headers
extern template void swap_if_needed(std::size_t& val, bool is_little_endian);
}  // namespace binary

inline std::string to_string(I32 x) {
#ifdef __INTEL_COMPILER
  /* Intel compilers with low GCC compat. don't have 32-bit specializations
   * for std::to_string, which makes them non-compliant with C++11
   */
  return std::to_string(static_cast<long long>(x));
#else
  return std::to_string(x);
#endif
}

}

#endif
