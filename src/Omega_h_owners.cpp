#include "Omega_h_owners.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_globals.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_linpart.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

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
    auto f = OMEGA_H_LAMBDA(LO old_owner) {
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
    parallel_for(nold_owners, f, "update_ownership(ranks)");
    copies2own_ranks = own_ranks;
  } else {
    Write<I32> old_owners2own_ranks(nold_owners);
    auto clients2ncopies = client2serv_comm->allgather(ncopies);
    auto f = OMEGA_H_LAMBDA(LO old_owner) {
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
    parallel_for(nold_owners, f, "update_ownership");
    copies2own_ranks =
        old_owners2copies.exch(Read<I32>(old_owners2own_ranks), 1);
  }
  auto copies2own_idxs =
      old_owners2copies.exch(Read<LO>(old_owners2own_idxs), 1);
  return Remotes(copies2own_ranks, copies2own_idxs);
}

/* uses update_ownership() with a linear partitioning of
   global numbers as the old partitioning.
   if global numbers obey good linear arrangement properties,
   the old partition should be decent and so this is an
   effective fallback if only globals are available */

Remotes owners_from_globals(
    CommPtr comm, Read<GO> globals, Read<I32> own_ranks) {
  auto copies2lins_dist = copies_to_linear_owners(comm, globals);
  return update_ownership(copies2lins_dist, own_ranks);
}

template <typename T>
Read<T> reduce_data_to_owners(
    Read<T> copy_data, Dist copies2owners, Int ncomps) {
  auto owners2copies = copies2owners.invert();
  auto serv_copy_data = copies2owners.exch(copy_data, ncomps);
  auto nowners = owners2copies.nroots();
  auto comm = copies2owners.parent_comm();
  auto owners2serv_copies = owners2copies.roots2items();
  auto owner_data_w = Write<T>(nowners * ncomps);
  auto f = OMEGA_H_LAMBDA(LO owner) {
    auto sc_begin = owners2serv_copies[owner];
    for (Int c = 0; c < ncomps; ++c) {
      owner_data_w[owner * ncomps + c] = serv_copy_data[sc_begin * ncomps + c];
    }
  };
  parallel_for(nowners, f, "reduce_data_to_owners");
  return owner_data_w;
}

GOs globals_from_owners(Mesh* mesh, Int ent_dim) {
  auto nnew_ents = mesh->nents(ent_dim);
  if (!mesh->could_be_shared(ent_dim)) {
    auto start = mesh->comm()->exscan(GO(nnew_ents), OMEGA_H_SUM);
    return Read<GO>(nnew_ents, start, 1);
  }
  auto new_owned = mesh->owned(ent_dim);
  auto new_globals = rescan_globals(mesh, new_owned);
  return mesh->sync_array(ent_dim, new_globals, 1);
}

#define INST(T)                                                                \
  template Read<T> reduce_data_to_owners(                                      \
      Read<T> copy_data, Dist copies2owners, Int ncomps);
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
