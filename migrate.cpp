Remotes form_old_use_owners(Mesh& old_mesh, Int high_dim, Int low_dim) {
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

LOs form_new_conn(Dist new_ents2old_owners, Dist old_owners2new_uses) {
  auto nnew_ents = new_ents2old_owners.nitems();
  auto serv_ents2new_idxs = new_ents2old_owners.exch(
      LOs(nnew_ents, 0, 1), 1);
  auto old_owners2new_ents = new_ents2old_owners.invert();
  auto serv_ents2ranks = old_owners2new_ents.items2ranks();
  auto serv_uses2ranks = old_owners2new_uses.items2ranks();
  auto nserv_uses = old_owners2new_uses.nitems();
  auto old_owners2serv_uses = old_owners2new_uses.roots2items();
  auto old_owners2serv_ents = old_owners2new_ents.roots2items();
  auto nold_owners = old_owners2new_ents.nroots();
  Write<LO> serv_uses2new_idxs(nserv_uses);
  auto f = LAMBDA(LO old_owner) {
    auto ebegin = old_owners2serv_ents[old_owner];
    auto eend = old_owners2serv_ents[old_owner + 1];
    auto ubegin = old_owners2serv_uses[old_owner];
    auto uend = old_owners2serv_uses[old_owner + 1];
    for (auto u = ubegin; u < uend; ++u) {
      auto rank = serv_uses2ranks[u];
      LO idx = -1;
      for (auto e = ebegin; e < eend; ++e) {
        if (serv_ents2ranks[e] == rank) {
          idx = serv_ents2new_idxs[e];
          break;
        }
      }
      serv_uses2new_idxs[u] = idx;
    }
  };
  parallel_for(nold_owners, f);
  auto serv_uses2new_uses = old_owners2new_uses;
  serv_uses2new_uses.set_roots2items(LOs());
  return serv_uses2new_uses.exch(LOs(serv_uses2new_idxs), 1);
}

void pull_down(Mesh& old_mesh, Int ent_dim, Int low_dim,
    Dist old_owners2new_ents,
    Adj& new_ents2new_lows, Dist& old_low_owners2new_lows) {
  auto nlows_per_high = simplex_degrees[ent_dim][low_dim];
  auto old_use_owners = form_old_use_owners(old_mesh,
      ent_dim, low_dim);
  auto new_use_own_ranks = old_owners2new_ents.exch(
      old_use_owners.ranks, nlows_per_high);
  auto new_use_own_idxs = old_owners2new_ents.exch(
      old_use_owners.idxs, nlows_per_high);
  Remotes new_use_owners(new_use_own_ranks, new_use_own_idxs);
  Dist low_uses2old_owners(old_mesh.comm(), new_use_owners,
      old_mesh.nents(low_dim));
  old_low_owners2new_lows = find_unique_use_owners(
      low_uses2old_owners);
  auto new_lows2old_owners = old_low_owners2new_lows.invert();
  auto old_low_owners2new_uses = low_uses2old_owners.invert();
  auto new_conn = form_new_conn(new_lows2old_owners, old_low_owners2new_uses);
  new_ents2new_lows.ab2b = new_conn;
  auto old_codes = old_mesh.ask_down(ent_dim, low_dim).codes;
  if (!old_codes.exists())
    return;
  auto new_codes = old_owners2new_ents.exch(old_codes, nlows_per_high);
  new_ents2new_lows.codes = new_codes;
}

void push_tags(Mesh const& old_mesh, Mesh& new_mesh,
    Int ent_dim, Dist old_owners2new_ents) {
  for (Int i = 0; i < old_mesh.ntags(ent_dim); ++i) {
    auto tag = old_mesh.get_tag(ent_dim, i);
    if (is<I8>(tag)) {
      auto array = to<I8>(tag)->array();
      array = old_owners2new_ents.exch(array, tag->ncomps());
      new_mesh.add_tag<I8>(ent_dim, tag->name(), tag->ncomps(), array);
    } else if (is<I32>(tag)) {
      auto array = to<I32>(tag)->array();
      array = old_owners2new_ents.exch(array, tag->ncomps());
      new_mesh.add_tag<I32>(ent_dim, tag->name(), tag->ncomps(), array);
    } else if (is<I64>(tag)) {
      auto array = to<I64>(tag)->array();
      array = old_owners2new_ents.exch(array, tag->ncomps());
      new_mesh.add_tag<I64>(ent_dim, tag->name(), tag->ncomps(), array);
    } else if (is<Real>(tag)) {
      auto array = to<Real>(tag)->array();
      array = old_owners2new_ents.exch(array, tag->ncomps());
      new_mesh.add_tag<Real>(ent_dim, tag->name(), tag->ncomps(), array);
    }
  }
}

void migrate_mesh(Mesh& old_mesh, Mesh& new_mesh,
    Remotes new_elems2old_owners) {
  auto comm = old_mesh.comm();
  auto dim = old_mesh.dim();
  new_mesh.set_comm(comm);
  new_mesh.set_dim(dim);
  auto nold_ents = old_mesh.nents(dim);
  Dist new_ents2old_owners(comm, new_elems2old_owners, nold_ents);
  auto old_owners2new_ents = new_ents2old_owners.invert();
  for (Int d = dim; d > VERT; --d) {
    Adj high2low;
    Dist old_low_owners2new_lows;
    pull_down(old_mesh, d, d - 1, old_owners2new_ents,
        high2low, old_low_owners2new_lows);
    new_mesh.set_ents(d, high2low);
    push_tags(old_mesh, new_mesh, dim, old_owners2new_ents);
    new_ents2old_owners = old_owners2new_ents.invert();
    auto owners = update_ownership(new_ents2old_owners, Read<I32>());
    new_mesh.set_tag<I32>(d, "owner", owners.ranks);
    new_mesh.set_own_idxs(d, owners.idxs);
    old_owners2new_ents = old_low_owners2new_lows;
  }
  auto new_verts2old_owners = old_owners2new_ents.invert();
  auto nnew_verts = new_verts2old_owners.nitems();
  new_mesh.set_verts(nnew_verts);
  push_tags(old_mesh, new_mesh, VERT, old_owners2new_ents);
  auto owners = update_ownership(new_ents2old_owners, Read<I32>());
  new_mesh.set_tag<I32>(VERT, "owner", owners.ranks);
  new_mesh.set_own_idxs(VERT, owners.idxs);
}

void migrate_mesh(Mesh& mesh, Remotes new_elems2old_owners) {
  Mesh new_mesh;
  migrate_mesh(mesh, new_mesh, new_elems2old_owners);
  mesh = new_mesh;
}
