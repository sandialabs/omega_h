void modify_conn(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    LOs prod_verts2verts,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs old_lows2new_lows) {
  auto down_degree = simplex_degrees[ent_dim][ent_dim - 1];
  auto old_ents2old_lows = old_mesh.ask_down(ent_dim, ent_dim - 1);
  auto old_ent_lows2old_lows = old_ents2old_lows.ab2b;
  auto old_ent_low_codes = old_ents2old_lows.codes;
  auto same_ent_lows2old_lows = unmap(
      same_ents2old_ents, old_ent_lows2old_lows, down_degree);
  auto same_ent_low_codes = unmap(
      same_ents2old_ents, old_ent_low_codes, 1);
  auto same_ent_lows2new_lows = compound_maps(
      same_ent_lows2old_lows, old_lows2new_lows);
  auto new_low_verts2new_verts = new_mesh.ask_verts_of(ent_dim - 1);
  auto new_verts2new_lows = new_mesh.ask_up(VERT, ent_dim - 1);
  auto prods2new_lows = reflect_down(prod_verts2verts,
      new_low_verts2new_verts, new_verts2new_lows, ent_dim, ent_dim - 1);
  auto prod_lows2new_lows = prods2new_lows.ab2b;
  auto prod_low_codes = prods2new_lows.codes;
  auto nsame_ents = same_ents2old_ents.size();
  auto nprods = prods2new_ents.size();
  auto nnew_ents = nsame_ents + nprods;
  Write<LO> new_ent_lows2new_lows(nnew_ents * down_degree);
  Write<I8> new_ent_low_codes(nnew_ents * down_degree);
  map_into(prod_lows2new_lows, prods2new_ents, new_ent_lows2new_lows,
      down_degree);
  map_into(same_ent_lows2new_lows, same_ents2new_ents, new_ent_lows2new_lows,
      down_degree);
  map_into(prod_low_codes, prods2new_ents, new_ent_low_codes,
      down_degree);
  map_into(same_ent_low_codes, same_ents2new_ents, new_ent_low_codes,
      down_degree);
  auto new_ents2new_lows = Adj(
      LOs(new_ent_lows2new_lows), Read<I8>(new_ent_low_codes));
  new_mesh.set_ents(ent_dim, new_ents2new_lows);
}

void modify_owners(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs old_ents2new_ents) {
  auto same_owners = unmap_owners(old_mesh, ent_dim,
      same_ents2old_ents, old_ents2new_ents);
  auto same_own_ranks = same_owners.ranks;
  auto same_own_idxs = same_owners.idxs;
  auto nprods = prods2new_ents.size();
  auto prod_own_ranks = Read<I32>(nprods, new_mesh.comm()->rank());
  auto prod_own_idxs = prods2new_ents;
  auto nsame_ents = same_ents2old_ents.size();
  auto nnew_ents = nsame_ents + nprods;
  Write<I32> new_own_ranks(nnew_ents);
  Write<LO> new_own_idxs(nnew_ents);
  map_into(same_own_ranks, same_ents2new_ents, new_own_ranks, 1);
  map_into(same_own_idxs, same_ents2new_ents, new_own_idxs, 1);
  map_into(prod_own_ranks, prods2new_ents, new_own_ranks, 1);
  map_into(prod_own_idxs, prods2new_ents, new_own_idxs, 1);
  auto new_owners = Remotes(Read<I32>(new_own_ranks), LOs(new_own_idxs));
  new_mesh.set_owners(ent_dim, new_owners);
}

LOs collect_same(Mesh& mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds) {
  auto nkds = mesh.nents(key_dim);
  auto kds_are_keys = mark_image(keys2kds, nkds);
  auto kds2ents = mesh.ask_graph(key_dim, ent_dim);
  auto kds2kd_ents = kds2ents.a2ab;
  auto kd_ents2ents = kds2ents.ab2b;
  auto kd_ents_are_adj = expand(kds_are_keys, kds2kd_ents, 1);
  auto ents_are_adj = permute(kd_ents_are_adj, kd_ents2ents, 1);
  auto ents_are_same = invert_marks(ents_are_adj);
  auto same_ents2old_ents = collect_marked(ents_are_same);
  return same_ents2old_ents;
}

LOs get_keys2reps(Mesh& mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds,
    LOs keys2nprods) {
  CHECK(ent_dim >= key_dim || ent_dim == VERT);
  auto nkeys = keys2kds.size();
  CHECK(nkeys == keys2nprods.size());
  LOs keys2reps;
  if (key_dim == ent_dim) {
    keys2reps = keys2kds;
  } else if (ent_dim > key_dim) {
    auto keys2ents2 = mesh.ask_up(key_dim, ent_dim);
    auto keys2key_ents = keys2ents2.a2ab;
    auto key_ents2ents = keys2ents2.ab2b;
    Write<LO> keys2reps_w;
    auto setup_reps = OSH_LAMBDA(LO key) {
      /* the first upward adjacent entity will represent
         this cavity during the updating of global numbers.
         upward ordering should be sorted by old globals
         already, so this is actually the adjacent entity
         with the lowest old global number */
      auto key_ent = keys2key_ents[key];
      auto rep = key_ents2ents[key_ent];
      keys2reps_w[key] = rep;
    };
    parallel_for(nkeys, setup_reps);
    keys2reps = keys2reps_w;
  } else {
    CHECK(ent_dim == VERT);
    CHECK(key_dim == EDGE);
    /* in the case of finding new globals for vertices after
       refining, we will use an endpoint vertex of the edge
       as the "representative" entity.
       this causes some concern because unlike before,
       when the reprentative was on the interior of the old cavity,
       it is now on the boundary.
       this is the reason for the atomic_add in get_local_rep_counts()
       and the exch_sum later on.
       I can't think of a nicer way to determine new vertex globals
       which is independent of partitioning and ordering */
    auto edge_verts2verts = mesh.ask_verts_of(EDGE);
    Write<LO> keys2reps_w;
    auto setup_reps = OSH_LAMBDA(LO key) {
      auto edge = keys2kds[key];
      keys2reps_w[key] = edge_verts2verts[edge * 2 + 0];
    };
    parallel_for(nkeys, setup_reps);
    keys2reps = keys2reps_w;
  }
  return keys2reps;
}

LOs get_rep_counts(
    Mesh& mesh,
    Int ent_dim,
    LOs keys2reps,
    LOs keys2nprods,
    LOs same_ents2ents) {
  auto nkeys = keys2reps.size();
  auto nents = mesh.nents(ent_dim);
  CHECK(nkeys == keys2nprods.size());
  auto nsame_ents = same_ents2ents.size();
  auto owned = mesh.owned(ent_dim);
  Write<LO> rep_counts(nents, 0);
  auto mark_same = OSH_LAMBDA(LO same_ent) {
    auto ent = same_ents2ents[same_ent];
    if (owned[ent])
      rep_counts[ent] = 1;
  };
  parallel_for(nsame_ents, mark_same);
  auto mark_reps = OSH_LAMBDA(LO key) {
    auto rep = keys2reps[key];
    auto nkey_prods = keys2nprods[key];
    atomic_add(&rep_counts[rep], nkey_prods);
  };
  parallel_for(nkeys, mark_reps);
  return rep_counts;
}

template <typename T>
void find_new_offsets(
    Int ent_dim,
    Read<T> old_ents2new_offsets,
    LOs same_ents2old_ents,
    LOs keys2reps,
    LOs keys2prods,
    Read<T>& same_ents2new_offsets,
    Read<T>& prods2new_offsets) {
  same_ents2new_offsets = unmap(same_ents2old_ents,
      old_ents2new_offsets, 1);
  auto keys2new_offsets = unmap(keys2reps,
      old_ents2new_offsets, 1);
  auto nprods = keys2prods.last();
  Write<T> prods2new_offsets_w(nprods);
  auto nkeys = keys2reps.size();
  CHECK(nkeys == keys2prods.size() - 1);
  auto write_prod_offsets = OSH_LAMBDA(LO key) {
    auto offset = keys2new_offsets[key];
    if (ent_dim == VERT) ++offset;
    for (auto prod = keys2prods[key]; prod < keys2prods[key]; ++prod) {
      prods2new_offsets_w[prod] = offset++;
    }
  };
  parallel_for(nkeys, write_prod_offsets);
  prods2new_offsets = prods2new_offsets_w;
}

template
void find_new_offsets(
    Int ent_dim,
    Read<LO> old_ents2new_offsets,
    LOs same_ents2old_ents,
    LOs keys2reps,
    LOs keys2prods,
    Read<LO>& same_ents2new_offsets,
    Read<LO>& prods2new_offsets);
template
void find_new_offsets(
    Int ent_dim,
    Read<GO> old_ents2new_offsets,
    LOs same_ents2old_ents,
    LOs keys2reps,
    LOs keys2prods,
    Read<GO>& same_ents2new_offsets,
    Read<GO>& prods2new_offsets);

void modify_globals(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs keys2reps,
    LOs rep_counts) {
  CHECK(ent_dim >= key_dim || ent_dim == VERT);
  auto nsame_ents = same_ents2old_ents.size();
  CHECK(nsame_ents == same_ents2new_ents.size());
  auto nkeys = keys2kds.size();
  CHECK(nkeys + 1 == keys2prods.size());
  auto nprods = prods2new_ents.size();
  CHECK(nprods == keys2prods.last());
  auto old_globals = old_mesh.get_array<GO>(ent_dim, "global");
  auto comm = old_mesh.comm();
  auto old_ents2lins = copies_to_linear_owners(comm, old_globals);
  auto lins2old_ents = old_ents2lins.invert();
  auto nlins = lins2old_ents.nitems();
  auto lin_rep_counts = old_ents2lins.exch_sum(LOs(rep_counts), 1);
  auto lin_local_offsets = offset_scan(lin_rep_counts);
  auto lin_global_count = lin_local_offsets.last();
  auto lin_global_offset = comm->exscan(lin_global_count, SUM);
  Write<GO> lin_globals(nlins);
  auto write_lin_globals = OSH_LAMBDA(LO lin) {
    lin_globals[lin] = lin_local_offsets[lin] + lin_global_offset;
  };
  parallel_for(nlins, write_lin_globals);
  auto old_ents2new_globals = lins2old_ents.exch(
      Read<GO>(lin_globals), 1);
  Read<GO> same_ents2new_globals;
  Read<GO> prods2new_globals;
  find_new_offsets(ent_dim, old_ents2new_globals, same_ents2old_ents,
      keys2reps, keys2prods, same_ents2new_globals, prods2new_globals);
  auto nnew_ents = new_mesh.nents(ent_dim);
  CHECK(nnew_ents == nsame_ents + nprods);
  Write<GO> new_globals(nnew_ents);
  map_into(same_ents2new_globals, same_ents2new_ents, new_globals, 1);
  map_into(prods2new_globals, prods2new_ents, new_globals, 1);
  new_mesh.add_tag(ent_dim, "global", 1, OSH_GLOBAL, Read<GO>(new_globals));
}
