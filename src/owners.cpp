#include "owners.hpp"

#include "linpart.hpp"
#include "loop.hpp"
#include "scan.hpp"

namespace osh {

Remotes update_ownership(Dist copies2old_owners, Read<I32> own_ranks) {
  auto ncopies = copies2old_owners.nitems();
  auto old_owners2copies = copies2old_owners.invert();
  auto nold_owners = old_owners2copies.nroots();
  auto serv_copies2copy_idxs = copies2old_owners.exch(LOs(ncopies, 0, 1), 1);
  auto client2serv_comm = copies2old_owners.comm();
  auto serv_copies2clients = old_owners2copies.items2msgs();
  auto old_owners2serv_copies = old_owners2copies.roots2items();
  auto clients2ranks = old_owners2copies.msgs2ranks();
  Write<LO> old_owners2own_idxs(nold_owners);
  Read<LO> copies2own_ranks;
  if (own_ranks.exists()) {
    auto serv_copies2own_ranks = copies2old_owners.exch(own_ranks, 1);
    auto f = LAMBDA(LO old_owner) {
      auto own_idx = -1;
      for (auto serv_copy = old_owners2serv_copies[old_owner];
           serv_copy < old_owners2serv_copies[old_owner + 1]; ++serv_copy) {
        auto client = serv_copies2clients[serv_copy];
        auto client_rank = clients2ranks[client];
        auto own_rank = serv_copies2own_ranks[serv_copy];
        if (own_rank == client_rank) {
          own_idx = serv_copies2copy_idxs[serv_copy];
          break;
        }
      }
      old_owners2own_idxs[old_owner] = own_idx;
    };
    parallel_for(nold_owners, f);
    copies2own_ranks = own_ranks;
  } else {
    Write<I32> old_owners2own_ranks(nold_owners);
    auto clients2ncopies = client2serv_comm->allgather(ncopies);
    auto f = LAMBDA(LO old_owner) {
      I32 own_rank = -1;
      LO nown_client_copies = -1;
      LO own_idx = -1;
      for (auto serv_copy = old_owners2serv_copies[old_owner];
           serv_copy < old_owners2serv_copies[old_owner + 1]; ++serv_copy) {
        auto client = serv_copies2clients[serv_copy];
        auto nclient_copies = clients2ncopies[client];
        auto client_rank = clients2ranks[client];
        if ((own_rank == -1) || (nclient_copies < nown_client_copies) ||
            ((nclient_copies == nown_client_copies) &&
             (client_rank < own_rank))) {
          auto copy_idx = serv_copies2copy_idxs[serv_copy];
          own_rank = client_rank;
          nown_client_copies = nclient_copies;
          own_idx = copy_idx;
        }
      }
      old_owners2own_ranks[old_owner] = own_rank;
      old_owners2own_idxs[old_owner] = own_idx;
    };
    parallel_for(nold_owners, f);
    copies2own_ranks =
        old_owners2copies.exch(Read<I32>(old_owners2own_ranks), 1);
  }
  auto copies2own_idxs =
      old_owners2copies.exch(Read<LO>(old_owners2own_idxs), 1);
  return Remotes(copies2own_ranks, copies2own_idxs);
}

Remotes owners_from_globals(CommPtr comm, Read<GO> globals,
                            Read<I32> own_ranks) {
  auto copies2lins_dist = copies_to_linear_owners(comm, globals);
  return update_ownership(copies2lins_dist, own_ranks);
}

template <typename T>
Read<T> reduce_data_to_owners(Read<T> copy_data, Dist copies2owners,
                              Int ncomps) {
  auto owners2copies = copies2owners.invert();
  auto serv_copy_data = copies2owners.exch(copy_data, ncomps);
  auto nowners = owners2copies.nroots();
  auto comm = copies2owners.parent_comm();
  auto owners2serv_copies = owners2copies.roots2items();
  auto owner_data_w = Write<T>(nowners * ncomps);
  auto f = LAMBDA(LO owner) {
    auto sc_begin = owners2serv_copies[owner];
    for (Int c = 0; c < ncomps; ++c) {
      owner_data_w[owner * ncomps + c] = serv_copy_data[sc_begin * ncomps + c];
    }
  };
  parallel_for(nowners, f);
  return owner_data_w;
}

void globals_from_owners(Mesh* new_mesh, Int ent_dim) {
  auto nnew_ents = new_mesh->nents(ent_dim);
  if (!new_mesh->could_be_shared(ent_dim)) {
    auto start = new_mesh->comm()->exscan(GO(nnew_ents), OSH_SUM);
    auto globals = Read<GO>(nnew_ents, start, 1);
    new_mesh->add_tag(ent_dim, "global", 1, OSH_GLOBAL, globals);
    return;
  }
  auto new_owned = new_mesh->owned(ent_dim);
  auto local_offsets = offset_scan(new_owned);
  auto nnew_owned = local_offsets.last();
  auto start = new_mesh->comm()->exscan(GO(nnew_owned), OSH_SUM);
  auto new_globals_w = Write<GO>(nnew_ents);
  parallel_for(nnew_ents,
               LAMBDA(LO e) { new_globals_w[e] = local_offsets[e] + start; });
  auto new_globals = Read<GO>(new_globals_w);
  new_globals = new_mesh->sync_array(ent_dim, new_globals, 1);
  new_mesh->add_tag(ent_dim, "global", 1, OSH_GLOBAL, new_globals);
}

#define INST(T)                                             \
  template Read<T> reduce_data_to_owners(Read<T> copy_data, \
                                         Dist copies2owners, Int ncomps);
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace osh
