#include "Omega_h_for.hpp"
#include "Omega_h_globals.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_mesh.hpp"

#include <iostream>

namespace Omega_h {

template <typename T>
GOs rescan_globals(Mesh* mesh, Read<T> counts, GO expected_total) {
  OMEGA_H_TIME_FUNCTION;
  auto const local_offsets = offset_scan(counts);
  auto const nnew = local_offsets.last();
  auto const start = mesh->comm()->exscan(GO(nnew), OMEGA_H_SUM);
  auto const new_globals_w = Write<GO>(counts.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    new_globals_w[i] = local_offsets[i] + start;
    if (expected_total != -1) {
      if (new_globals_w[i] > expected_total) {
        std::cerr << "BAD new_globals_w " << new_globals_w[i] << " local_offsets[" << i << "] "
          << local_offsets[i] << " start " << start << " expected total " << expected_total << '\n';
      }
    //OMEGA_H_CHECK(new_globals_w[i] <= expected_total);
    }
  };
  parallel_for(counts.size(), std::move(f));
  return new_globals_w;
}

#define OMEGA_H_INST_DECL(T) \
template GOs rescan_globals(Mesh*, Read<T>, GO);
OMEGA_H_INST_DECL(I8)
OMEGA_H_INST_DECL(I32)
#undef OMEGA_H_INST_DECL

}  // end namespace Omega_h
