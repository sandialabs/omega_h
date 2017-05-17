#include "Omega_h_bipart.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_map.hpp"

namespace Omega_h {

Dist bi_partition(CommPtr comm, Read<I8> marks) {
  auto halfsize = divide_no_remainder(comm->size(), 2);
  Write<I32> dest_ranks(marks.size());
  Write<LO> dest_idxs(marks.size());
  LO linsize = -1;
  I32 rank_start = 0;
  for (Int half = 0; half < 2; ++half) {
    marks = invert_marks(marks);
    auto marked = collect_marked(marks);
    auto total = comm->allreduce(GO(marked.size()), OMEGA_H_SUM);
    auto start = comm->exscan(GO(marked.size()), OMEGA_H_SUM);
    Read<GO> globals(marked.size(), start, 1);
    auto owners = globals_to_linear_owners(globals, total, halfsize);
    map_into(add_to_each(owners.ranks, rank_start), marked, dest_ranks, 1);
    map_into(owners.idxs, marked, dest_idxs, 1);
    if (rank_start <= comm->rank() && comm->rank() < (rank_start + halfsize)) {
      linsize =
          linear_partition_size(total, halfsize, comm->rank() - rank_start);
    }
    rank_start += halfsize;
  }
  auto dests = Remotes(Read<I32>(dest_ranks), Read<LO>(dest_idxs));
  return Dist(comm, dests, linsize);
}

}  // end namespace Omega_h
