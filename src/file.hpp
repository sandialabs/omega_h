#ifndef FILE_HPP
#define FILE_HPP

#include "internal.hpp"

namespace osh {

bool is_little_endian_cpu();
void safe_mkdir(const char* path);
bool directory_exists(const char* path);
std::string parent_path(std::string const& path);
std::string path_leaf_name(std::string const& path);

namespace binary {

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
void read(std::istream& stream, Mesh* mesh);

#define INST_DECL(T)                                                     \
  extern template void swap_if_needed(T& val, bool is_little_endian);    \
  extern template Read<T> swap_if_needed(Read<T> array,                  \
                                         bool is_little_endian);         \
  extern template void write_value(std::ostream& stream, T val);         \
  extern template void read_value(std::istream& stream, T& val);         \
  extern template void write_array(std::ostream& stream, Read<T> array); \
  extern template void read_array(std::istream& stream, Read<T>& array,  \
                                  bool is_compressed);
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL
// for VTK compression headers
extern template void swap_if_needed(std::size_t& val, bool is_little_endian);
}

}  // end namespace osh

#endif
