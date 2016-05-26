Remotes owners_from_globals(CommPtr comm,
    Read<GO> globals, Read<I32> own_ranks) {
  auto ncopies = globals.size();
  auto total = find_total_globals(comm, globals);
  auto nlins = linear_partition_size(comm, total);
  auto copies2lins_map = globals_to_linear_owners(comm, globals, total);
  auto copies2lins_dist = Dist(comm, copies2lins_map, globals.size());
  auto serv_copies2copy_idxs = copies2lins_dist.exch(LOs(ncopies, 0, 1), 1);
  auto client2serv_comm = copies2lins_dist.comm();
  auto lins2copies_dist = copies2lins_dist.invert();
  auto serv_copies2clients = lins2copies_dist.items2msgs();
  auto lins2serv_copies = lins2copies_dist.roots2items();
  auto clients2ranks = lins2copies_dist.msgs2ranks();
  Write<LO> lins2own_idxs(nlins);
  Read<LO> copies2own_ranks;
  if (own_ranks.exists()) {
    auto serv_copies2own_ranks = copies2lins_dist.exch(own_ranks, 1);
    auto f = LAMBDA(LO lin) {
      for (auto serv_copy = lins2serv_copies[lin];
           serv_copy < lins2serv_copies[lin + 1];
           ++serv_copy) {
        auto client = serv_copies2clients[serv_copy];
        auto client_rank = clients2ranks[client];
        auto own_rank = serv_copies2own_ranks[serv_copy];
        if (own_rank == client_rank) {
          auto own_idx = serv_copies2copy_idxs[serv_copy];
          lins2own_idxs[lin] = own_idx;
          break;
        }
      }
    };
    parallel_for(nlins, f);
    copies2own_ranks = own_ranks;
  } else {
    Write<I32> lins2own_ranks(nlins);
    auto clients2ncopies = client2serv_comm->allgather(ncopies);
    auto f = LAMBDA(LO lin) {
      I32 own_rank = -1;
      LO nown_client_copies = -1;
      LO own_idx = -1;
      for (auto serv_copy = lins2serv_copies[lin];
           serv_copy < lins2serv_copies[lin + 1];
           ++serv_copy) {
        auto client = serv_copies2clients[serv_copy];
        auto nclient_copies = clients2ncopies[client];
        auto client_rank = clients2ranks[client];
        if ((own_rank == -1) ||
            (nclient_copies < nown_client_copies) ||
            ((nclient_copies == nown_client_copies) &&
             (client_rank < own_rank))) {
          auto copy_idx = serv_copies2copy_idxs[serv_copy];
          own_rank = client_rank;
          nown_client_copies = nclient_copies;
          own_idx = copy_idx;
        }
      }
      lins2own_ranks[lin] = own_rank;
      lins2own_idxs[lin] = own_idx;
    };
    parallel_for(nlins, f);
    copies2own_ranks = lins2copies_dist.exch(Read<I32>(lins2own_ranks), 1);
  }
  auto copies2own_idxs = lins2copies_dist.exch(Read<LO>(lins2own_idxs), 1);
  return Remotes(copies2own_ranks, copies2own_idxs);
}
