#include "Omega_h_for.hpp"
#include "Omega_h_globals.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

template <typename T>
GOs rescan_globals(Mesh* mesh, Read<T> counts) {
  begin_code("rescan_globals");
  auto local_offsets = offset_scan(counts);
  auto nnew = local_offsets.last();
  auto start = mesh->comm()->exscan(GO(nnew), OMEGA_H_SUM);
  auto new_globals_w = Write<GO>(counts.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    new_globals_w[i] = local_offsets[i] + start;
  };
  parallel_for(counts.size(), f, "rescan_globals");
  return new_globals_w;
}

#define OMEGA_H_INST_DECL(T) \
template GOs rescan_globals(Mesh*, Read<T>);
OMEGA_H_INST_DECL(I8)
OMEGA_H_INST_DECL(I32)
#undef OMEGA_H_INST_DECL

}  // end namespace Omega_h
