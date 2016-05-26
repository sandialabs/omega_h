Remotes form_down_use_owners(Mesh& old_mesh, Int high_dim, Int low_dim) {
  auto uses2lows = old_mesh.ask_down(high_dim, low_dim).ab2b;
  auto lows2owners = old_mesh.ask_owners(low_dim);
  auto own_ranks = unmap(uses2lows, lows2owners.ranks, 1);
  auto own_idxs = unmap(uses2lows, lows2owners.idxs, 1);
  return Remotes(own_ranks, own_idxs);
}

Dist find_unique_use_owners(Dist uses2old_owners) {
  auto old_owners2uses = uses2old_owners.invert();
  auto nold_owners = old_owners2uses.nroots();
  auto nserv_uses = old_owners2uses.nitems();
  auto serv_uses2ranks = old_owners2uses.items2ranks();
  auto old_owners2serv_uses = old_owners2uses.roots2items();
  Write<I8> keep(nserv_uses, 1);
  Write<LO> degrees(nold_owners);
  auto f = LAMBDA(LO old_owner) {
    LO degree = 0;
    auto begin = old_owners2serv_uses[old_owner];
    auto end = old_owners2serv_uses[old_owner + 1];
    for (LO serv_use = begin; serv_use < end; ++serv_use) {
      if (!keep[serv_use])
        continue;
      ++degree;
      auto rank = serv_uses2ranks[serv_use];
      for (LO serv_use2 = serv_use + 1; serv_use2 < end; ++serv_use2) {
        if (serv_uses2ranks[serv_use2] == rank) {
          keep[serv_use2] = 0;
        }
      }
    }
    degrees[old_owner] = degree;
  };
  parallel_for(nold_owners, f);
  auto uniq_serv_uses = collect_marked(Read<I8>(keep));
  auto uniq_serv_uses2ranks = unmap(uniq_serv_uses, serv_uses2ranks, 1);
  auto old_owners2uniq_serv_uses = offset_scan(LOs(degrees));
  Dist old_owners2uniq_uses;
  old_owners2uniq_uses.set_parent_comm(uses2old_owners.parent_comm());
  old_owners2uniq_uses.set_dest_ranks(uniq_serv_uses2ranks);
  old_owners2uniq_uses.set_roots2items(old_owners2uniq_serv_uses);
  return old_owners2uniq_uses;
}
