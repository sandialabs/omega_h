#include "linpart.hpp"

#include "array.hpp"
#include "loop.hpp"

namespace osh {

Remotes globals_to_linear_owners(Read<GO> globals, GO total, I32 comm_size) {
  auto comm_size_gt = GO(comm_size);
  auto quot = total / comm_size_gt;
  auto rem = total % comm_size_gt;
  auto split = ((quot + 1) * rem);
  Write<I32> ranks(globals.size());
  Write<LO> idxs(globals.size());
  auto f = LAMBDA(LO i) {
    if (globals[i] < split) {
      ranks[i] = static_cast<I32>(globals[i] / (quot + 1));
      idxs[i] = static_cast<LO>(globals[i] % (quot + 1));
    } else {
      auto g = globals[i] - split;
      ranks[i] = static_cast<I32>(g / quot + rem);
      idxs[i] = static_cast<LO>(g % quot);
    }
  };
  parallel_for(globals.size(), f);
  return Remotes(ranks, idxs);
}

LO linear_partition_size(GO total, I32 comm_size, I32 comm_rank) {
  auto comm_size_gt = GO(comm_size);
  auto quot = total / comm_size_gt;
  auto rem = total % comm_size_gt;
  if (comm_rank < static_cast<I32>(rem))
    return static_cast<LO>(quot + 1);
  else
    return static_cast<LO>(quot);
}

Remotes globals_to_linear_owners(CommPtr comm, Read<GO> globals, GO total) {
  return globals_to_linear_owners(globals, total, comm->size());
}

LO linear_partition_size(CommPtr comm, GO total) {
  return linear_partition_size(total, comm->size(), comm->rank());
}

GO find_total_globals(CommPtr comm, Read<GO> globals) {
  auto a = comm->allreduce(max(globals), OSH_MAX);
  if (a < 0) return 0;
  return a + 1;
}

Dist copies_to_linear_owners(CommPtr comm, Read<GO> globals) {
  auto total = find_total_globals(comm, globals);
  auto nlins = linear_partition_size(comm, total);
  auto copies2lins_map = globals_to_linear_owners(comm, globals, total);
  auto copies2lins_dist = Dist(comm, copies2lins_map, nlins);
  return copies2lins_dist;
}

}  // end namespace osh
