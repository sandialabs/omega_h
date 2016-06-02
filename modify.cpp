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

void modify_globals(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2ents,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  CHECK(ent_dim >= key_dim || ent_dim == VERT);
  auto nold_ents = old_mesh.nents(ent_dim);
  auto nsame_ents = same_ents2old_ents.size();
  CHECK(nsame_ents == same_ents2new_ents.size());
  auto nkeys = keys2ents.size();
  CHECK(nkeys + 1 == keys2prods.size());
  auto nprods = prods2new_ents.size();
  CHECK(nprods == keys2prods.last());
  auto old_owned = old_mesh.owned(ent_dim);
  Write<LO> old_locals(nold_ents, 0);
  auto mark_old_same = LAMBDA(LO same_ent) {
    auto old_ent = same_ents2old_ents[same_ent];
    if (old_owned[old_ent])
      old_locals[old_ent] = 1;
  };
  parallel_for(nsame_ents, mark_old_same);
  auto keys2nprods = get_degrees(keys2prods);
  LOs keys2reps;
  if (key_dim == ent_dim) {
    keys2reps = keys2ents;
  } else if (ent_dim > key_dim) {
    auto keys2ents2 = old_mesh.ask_up(key_dim, ent_dim);
    auto keys2key_ents = keys2ents2.a2ab;
    auto key_ents2ents = keys2ents2.ab2b;
    Write<LO> keys2reps_w;
    auto setup_reps = LAMBDA(LO key) {
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
       this is the reason for the atomic_add for old_locals
       and the exch_sum later on.
       I can't think of a nicer way to determine new vertex globals
       which is independent of partitioning and ordering */
    auto edge_verts2verts = old_mesh.ask_verts_of(EDGE);
    Write<LO> keys2reps_w;
    auto setup_reps = LAMBDA(LO key) {
      auto edge = keys2ents[key];
      keys2reps_w[key] = edge_verts2verts[edge * 2 + 0];
    };
    parallel_for(nkeys, setup_reps);
    keys2reps = keys2reps_w;
  }
  auto mark_reps = LAMBDA(LO key) {
    auto rep = keys2reps[key];
    auto nkey_prods = keys2nprods[key];
    atomic_add(&old_locals[rep], nkey_prods);
  };
  parallel_for(nkeys, mark_reps);
  auto old_globals = old_mesh.get_array<GO>(ent_dim, "global");
  auto comm = old_mesh.comm();
  auto old_ents2lins = copies_to_linear_owners(comm, old_globals);
  auto lins2old_ents = old_ents2lins.invert();
  auto nlins = lins2old_ents.nitems();
  auto lin_locals = old_ents2lins.exch_sum(LOs(old_locals), 1);
  auto lin_local_offsets = offset_scan(lin_locals);
  auto lin_global_count = lin_local_offsets.last();
  auto lin_global_offset = comm->exscan(lin_global_count, SUM);
  Write<GO> lin_globals(nlins);
  auto write_lin_globals = LAMBDA(LO lin) {
    lin_globals[lin] = lin_local_offsets[lin] + lin_global_offset;
  };
  parallel_for(nlins, write_lin_globals);
  auto old_ents2new_globals = lins2old_ents.exch(
      Read<GO>(lin_globals), 1);
  auto same_ents2new_globals = unmap(same_ents2old_ents,
      old_ents2new_globals, 1);
  auto keys2rep_globals = unmap(keys2reps, old_ents2new_globals, 1);
  auto nnew_ents = new_mesh.nents(ent_dim);
  CHECK(nnew_ents == nsame_ents + nprods);
  Write<GO> new_globals(nnew_ents);
  map_into(same_ents2new_globals, same_ents2new_ents, new_globals, 1);
  auto write_cavity_globals = LAMBDA(LO key) {
    auto global = keys2rep_globals[key];
    if (ent_dim == VERT) ++global;
    for (auto prod = keys2prods[key]; prod < keys2prods[key]; ++prod) {
      auto new_ent = prods2new_ents[prod];
      new_globals[new_ent] = global++;
    }
  };
  parallel_for(nkeys, write_cavity_globals);
  new_mesh.add_tag(ent_dim, "global", 1, Read<GO>(new_globals));
}
