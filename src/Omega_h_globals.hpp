#ifndef OMEGA_H_GLOBALS_HPP
#define OMEGA_H_GLOBALS_HPP

#include "Omega_h_array.hpp"

namespace Omega_h {

class Mesh;

template <typename T>
GOs rescan_globals(Mesh* mesh, Read<T> counts);

#define OMEGA_H_INST_DECL(T) extern template GOs rescan_globals(Mesh*, Read<T>);
OMEGA_H_INST_DECL(I8)
OMEGA_H_INST_DECL(I32)
#undef OMEGA_H_INST_DECL

}  // end namespace Omega_h

#endif
