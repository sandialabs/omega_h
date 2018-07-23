#ifndef OMEGA_H_FILE_HPP
#define OMEGA_H_FILE_HPP

#include <iosfwd>
#include <vector>

#include <Omega_h_config.h>
#include <Omega_h_array.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_tag.hpp>

namespace Omega_h {

bool ends_with(std::string const& s, std::string const& suffix);
bool is_little_endian_cpu();
void safe_mkdir(const char* path);
bool file_exists(const char* path);
bool directory_exists(const char* path);
std::string parent_path(std::string const& path);
std::string path_leaf_name(std::string const& path);

#ifdef OMEGA_H_USE_LIBMESHB
namespace meshb {
void read(Mesh* mesh, std::string const& filepath);
void write(Mesh* mesh, std::string const& filepath, int version = 2);
void read_sol(
    Mesh* mesh, std::string const& filepath, std::string const& sol_name);
void write_sol(Mesh* mesh, std::string const& filepath,
    std::string const& sol_name, int version = 2);
}  // namespace meshb
#endif

#ifdef OMEGA_H_USE_SEACASEXODUS
namespace exodus {
enum ClassifyWith {
  NODE_SETS = 0x1,
  SIDE_SETS = 0x2,
};
int open(std::string const& path, bool verbose = false);
void close(int exodus_file);
int get_num_time_steps(int exodus_file);
void read_mesh(int exodus_file, Mesh* mesh, bool verbose = false,
    int classify_with = NODE_SETS | SIDE_SETS);
void read_nodal_fields(int exodus_file, Mesh* mesh, int time_step,
    std::string const& prefix = "", std::string const& postfix = "",
    bool verbose = false);
void write(std::string const& path, Mesh* mesh, bool verbose = false,
    int classify_with = NODE_SETS | SIDE_SETS);
Mesh read_sliced(std::string const& path, CommPtr comm, bool verbose = false,
    int classify_with = NODE_SETS | SIDE_SETS, int time_step = -1);
}  // namespace exodus
#endif

namespace gmsh {
Mesh read(std::istream& stream, CommPtr comm);
Mesh read(std::string const& filename, CommPtr comm);
void write(std::ostream& stream, Mesh* mesh);
void write(std::string const& filepath, Mesh* mesh);
}  // namespace gmsh

namespace vtk {
static constexpr bool do_compress = true;
static constexpr bool dont_compress = false;
TagSet get_all_vtk_tags(Mesh* mesh, Int cell_dim);
void write_vtu(std::ostream& stream, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress = true);
void write_vtu(std::string const& filename, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress = true);
void write_vtu(std::string const& filename, Mesh* mesh, Int cell_dim,
    bool compress = true);
void write_vtu(std::string const& filename, Mesh* mesh, bool compress = true);
void write_parallel(std::string const& path, Mesh* mesh, Int cell_dim,
    TagSet const& tags, bool compress = true);
void write_parallel(
    std::string const& path, Mesh* mesh, Int cell_dim, bool compress = true);
void write_parallel(std::string const& path, Mesh* mesh, bool compress = true);
class Writer {
  Mesh* mesh_;
  std::string root_path_;
  Int cell_dim_;
  bool compress_;
  I64 step_;
  std::streampos pvd_pos_;

 public:
  Writer();
  Writer(Writer const&) = default;
  Writer& operator=(Writer const&) = default;
  ~Writer() = default;
  Writer(std::string const& root_path, Mesh* mesh, Int cell_dim = -1,
      Real restart_time = 0.0, bool compress = true);
  void write();
  void write(Real time);
  void write(Real time, TagSet const& tags);
  void write(I64 step, Real time, TagSet const& tags);
};
class FullWriter {
  std::vector<Writer> writers_;

 public:
  FullWriter() = default;
  FullWriter(std::string const& root_path, Mesh* mesh, Real restart_time = 0.0,
      bool compress = true);
  void write(Real time);
  void write();
};
}  // end namespace vtk

namespace binary {

void write(std::string const& path, Mesh* mesh);
Mesh read(std::string const& path, Library* lib, bool strict = false);
Mesh read(std::string const& path, CommPtr comm, bool strict = false);
I32 read(
    std::string const& path, CommPtr comm, Mesh* mesh, bool strict = false);
I32 read_nparts(std::string const& path, CommPtr comm);
I32 read_version(std::string const& path, CommPtr comm);
void read_in_comm(
    std::string const& path, CommPtr comm, Mesh* mesh, I32 version);

constexpr I32 latest_version = 8;

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
extern template void swap_if_needed(std::uint64_t& val, bool is_little_endian);
}  // namespace binary

void write_reals_txt(std::string const& filename, Reals a, Int ncomps);
void write_reals_txt(std::ostream& stream, Reals a, Int ncomps);
Reals read_reals_txt(std::string const& filename, LO n, Int ncomps);
Reals read_reals_txt(std::istream& stream, LO n, Int ncomps);

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

}  // namespace Omega_h

#endif
