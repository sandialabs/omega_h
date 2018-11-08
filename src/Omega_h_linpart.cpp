#include "Omega_h_linpart.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"

namespace Omega_h {

Remotes globals_to_linear_owners(Read<GO> globals, GO total, I32 comm_size) {
  auto comm_size_gt = GO(comm_size);
  auto quot = total / comm_size_gt;
  auto rem = total % comm_size_gt;
  auto split = ((quot + 1) * rem);
  Write<I32> ranks(globals.size());
  Write<LO> idxs(globals.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    if (globals[i] < split) {
      ranks[i] = static_cast<I32>(globals[i] / (quot + 1));
      idxs[i] = static_cast<LO>(globals[i] % (quot + 1));
    } else {
      auto g = globals[i] - split;
      ranks[i] = static_cast<I32>(g / quot + rem);
      idxs[i] = static_cast<LO>(g % quot);
    }
  };
  parallel_for(globals.size(), std::move(f));
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
  auto a = get_max(comm, globals);
  if (a < 0) return 0;
  return a + 1;
}

Dist copies_to_linear_owners(CommPtr comm, Read<GO> globals, GO expected_total) {
  OMEGA_H_TIME_FUNCTION;
  auto total = find_total_globals(comm, globals);
  if (expected_total != -1) OMEGA_H_CHECK(total == expected_total);
  auto nlins = linear_partition_size(comm, total);
  auto copies2lins_map = globals_to_linear_owners(comm, globals, total);
  if (total % GO(nlins) == 0) {
    for (int i = 0; i < copies2lins_map.idxs.size(); ++i) {
      OMEGA_H_CHECK(0 <= copies2lins_map.idxs[i]);
      OMEGA_H_CHECK(copies2lins_map.idxs[i] < nlins);
    }
  }
  auto copies2lins_dist = Dist(comm, copies2lins_map, nlins, (total % GO(nlins) == 0));
  return copies2lins_dist;
}

}  // end namespace Omega_h
