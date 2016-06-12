static void modify_conn(Mesh* old_mesh, Mesh* new_mesh,
    Int ent_dim,
    LOs prod_verts2verts,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs old_lows2new_lows) {
  auto low_dim = ent_dim - 1;
  auto down_degree = simplex_degrees[ent_dim][low_dim];
  auto old_ents2old_lows = old_mesh->ask_down(ent_dim, low_dim);
  auto old_ent_lows2old_lows = old_ents2old_lows.ab2b;
  auto same_ent_lows2old_lows = unmap(
      same_ents2old_ents, old_ent_lows2old_lows, down_degree);
  auto same_ent_lows2new_lows = compound_maps(
      same_ent_lows2old_lows, old_lows2new_lows);
  auto prods2new_lows = Adj();
  if (low_dim > VERT) {
    auto new_low_verts2new_verts = new_mesh->ask_verts_of(low_dim);
    auto new_verts2new_lows = new_mesh->ask_up(VERT, low_dim);
    prods2new_lows = reflect_down(prod_verts2verts,
      new_low_verts2new_verts, new_verts2new_lows, ent_dim, low_dim);
  } else {
    prods2new_lows = Adj(prod_verts2verts);
  }
  auto prod_lows2new_lows = prods2new_lows.ab2b;
  auto nsame_ents = same_ents2old_ents.size();
  auto nprods = prods2new_ents.size();
  auto nnew_ents = nsame_ents + nprods;
  Write<LO> new_ent_lows2new_lows(nnew_ents * down_degree);
  auto new_ent_low_codes = Write<I8>();
  map_into(prod_lows2new_lows, prods2new_ents, new_ent_lows2new_lows,
      down_degree);
  map_into(same_ent_lows2new_lows, same_ents2new_ents, new_ent_lows2new_lows,
      down_degree);
  if (low_dim > VERT) {
    new_ent_low_codes = Write<I8>(nnew_ents * down_degree);
    auto old_ent_low_codes = old_ents2old_lows.codes;
    CHECK(same_ents2old_ents.size() == same_ents2new_ents.size());
    auto same_ent_low_codes = unmap(
        same_ents2old_ents, old_ent_low_codes, down_degree);
    map_into(same_ent_low_codes, same_ents2new_ents, new_ent_low_codes,
        down_degree);
    auto prod_low_codes = prods2new_lows.codes;
    map_into(prod_low_codes, prods2new_ents, new_ent_low_codes,
        down_degree);
  }
  auto new_ents2new_lows = Adj(
      LOs(new_ent_lows2new_lows), Read<I8>(new_ent_low_codes));
  new_mesh->set_ents(ent_dim, new_ents2new_lows);
}

static void modify_owners(Mesh* old_mesh, Mesh* new_mesh,
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
  auto prod_own_ranks = Read<I32>(nprods, new_mesh->comm()->rank());
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
  new_mesh->set_owners(ent_dim, new_owners);
}

static LOs collect_same(Mesh* mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds) {
  if (ent_dim < key_dim) {
    auto nents = mesh->nents(ent_dim);
    return LOs(nents, 0, 1);
  }
  auto nkds = mesh->nents(key_dim);
  auto kds_are_keys = mark_image(keys2kds, nkds);
  auto ents_are_adj = Read<I8>();
  if (ent_dim == key_dim) {
    ents_are_adj = kds_are_keys;
  } else {
    CHECK(ent_dim > key_dim);
    ents_are_adj = mark_up(mesh, key_dim, ent_dim, kds_are_keys);
  }
  auto ents_not_adj = invert_marks(ents_are_adj);
  return collect_marked(ents_not_adj);
}

static LOs get_keys2reps(Mesh* mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds,
    LOs keys2nprods) {
  CHECK(ent_dim >= key_dim || (ent_dim == VERT && key_dim == EDGE));
  auto nkeys = keys2kds.size();
  CHECK(nkeys == keys2nprods.size());
  LOs keys2reps;
  if (key_dim == ent_dim) {
    keys2reps = keys2kds;
  } else if (ent_dim > key_dim) {
    auto kds2ents = mesh->ask_up(key_dim, ent_dim);
    auto kds2kd_ents = kds2ents.a2ab;
    auto kd_ents2ents = kds2ents.ab2b;
    Write<LO> keys2reps_w(nkeys);
    auto setup_reps = LAMBDA(LO key) {
      /* the first upward adjacent entity will represent
         this cavity during the updating of global numbers.
         upward ordering should be sorted by old globals
         already, so this is actually the adjacent entity
         with the lowest old global number */
      auto kd = keys2kds[key];
      auto first_kd_ent = kds2kd_ents[kd];
      auto first_adj = kd_ents2ents[first_kd_ent];
      auto rep = first_adj;
      keys2reps_w[key] = rep;
    };
    parallel_for(nkeys, setup_reps);
    keys2reps = keys2reps_w;
  } else {
    CHECK(ent_dim == VERT && key_dim == EDGE);
    /* in the case of finding new globals for vertices after
       refining, we will use an endpoint vertex of the edge
       as the "representative" entity.
       this causes some concern because unlike before,
       when the reprentative was on the interior of the old cavity,
       it is now on the boundary.
       this is the reason for the atomic_add in get_rep_counts()
       and the exch_sum later on.
       I can't think of a nicer way to determine new vertex globals
       which is independent of partitioning and ordering */
    auto edge_verts2verts = mesh->ask_verts_of(EDGE);
    Write<LO> keys2reps_w(nkeys);
    auto setup_reps = LAMBDA(LO key) {
      auto edge = keys2kds[key];
      keys2reps_w[key] = edge_verts2verts[edge * 2 + 0];
    };
    parallel_for(nkeys, setup_reps);
    keys2reps = keys2reps_w;
  }
  return keys2reps;
}

static LOs get_rep_counts(
    Mesh* mesh,
    Int ent_dim,
    LOs keys2reps,
    LOs keys2nprods,
    LOs same_ents2ents,
    bool are_global) {
  auto nkeys = keys2reps.size();
  auto nents = mesh->nents(ent_dim);
  CHECK(nkeys == keys2nprods.size());
  auto nsame_ents = same_ents2ents.size();
  auto owned = mesh->owned(ent_dim);
  CHECK(owned.size() == nents);
  Write<LO> rep_counts(nents, 0);
  auto mark_same = LAMBDA(LO same_ent) {
    auto ent = same_ents2ents[same_ent];
    CHECK(ent < nents);
    if ((!are_global) || owned[ent])
      rep_counts[ent] = 1;
  };
  parallel_for(nsame_ents, mark_same);
  auto mark_reps = LAMBDA(LO key) {
    auto rep = keys2reps[key];
    auto nkey_prods = keys2nprods[key];
    atomic_add(&rep_counts[rep], nkey_prods);
  };
  parallel_for(nkeys, mark_reps);
  return rep_counts;
}

/* one more thing we need to for the pathological
   case of new vertices is to resolve the situation
   where multiple cavities share the same representative
   vertex of the old mesh->
   that vertex represents the sum of all adjacent cavities
   and has a single global number from which the new globals
   of all new cavity vertices are computed.
   we just need to determine an order in which to assign
   adjacent globals from the representative.
   the most intuitive order is the upward adjacent ordering
   from vertices to edges (guaranteed to be sorted by globals).
   however, this is only available in full in GHOSTED mode.
   so, this function writes down that ordering while in GHOSTED
   mode so it can be used later by find_new_offsets, which
   runs in ELEMENT_BASED mode */

LOs get_edge2rep_order(Mesh* mesh, Read<I8> edges_are_keys) {
  auto nedges = mesh->nents(EDGE);
  auto nverts = mesh->nents(VERT);
  auto order_w = Write<LO>(nedges, -1);
  auto verts2edges = mesh->ask_up(VERT, EDGE);
  auto v2ve = verts2edges.a2ab;
  auto ve2e = verts2edges.ab2b;
  auto ve_codes = verts2edges.codes;
  auto f = LAMBDA(LO v) {
    LO i = 0;
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto code = ve_codes[ve];
      auto dir = code_which_down(code);
      if (dir == 0 && edges_are_keys[e]) {
        order_w[e] = i++;
      }
    }
  };
  parallel_for(nverts, f);
  return order_w;
}

template <typename T>
static void find_new_offsets(
    Read<T> old_ents2new_offsets,
    LOs same_ents2old_ents,
    LOs keys2kds,
    LOs keys2reps,
    LOs keys2prods,
    LOs edge2rep_order,
    Read<T>* p_same_ents2new_offsets,
    Read<T>* p_prods2new_offsets) {
  *p_same_ents2new_offsets = unmap(same_ents2old_ents,
      old_ents2new_offsets, 1);
  auto keys2new_offsets = unmap(keys2reps,
      old_ents2new_offsets, 1);
  auto nprods = keys2prods.last();
  Write<T> prods2new_offsets_w(nprods);
  auto nkeys = keys2reps.size();
  CHECK(nkeys == keys2prods.size() - 1);
  if (edge2rep_order.exists()) {
    CHECK(keys2kds.exists());
    auto write_prod_offsets = LAMBDA(LO key) {
      // plus one because the representatives
      // exist in the new mesh, so they will get
      // the global number of that vertex in the new mesh,
      // which gets saved to same_ents2new_offsets above
      // the globals for new vertices start after that one,
      // and are ordered by edge2rep_order
      auto offset = keys2new_offsets[key] + 1;
      auto edge = keys2kds[key];
      auto prod = key;
      prods2new_offsets_w[prod] = offset + edge2rep_order[edge];
    };
    parallel_for(nkeys, write_prod_offsets);
  } else {
    auto write_prod_offsets = LAMBDA(LO key) {
      auto offset = keys2new_offsets[key];
      for (auto prod = keys2prods[key]; prod < keys2prods[key + 1]; ++prod) {
        prods2new_offsets_w[prod] = offset;
        ++offset;
      }
    };
    parallel_for(nkeys, write_prod_offsets);
  }
  *p_prods2new_offsets = prods2new_offsets_w;
}

static void modify_globals(Mesh* old_mesh, Mesh* new_mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs keys2reps,
    LOs global_rep_counts) {
  CHECK(ent_dim >= key_dim || (ent_dim == VERT && key_dim == EDGE));
  auto nsame_ents = same_ents2old_ents.size();
  CHECK(nsame_ents == same_ents2new_ents.size());
  auto nkeys = keys2kds.size();
  CHECK(nkeys + 1 == keys2prods.size());
  auto nprods = prods2new_ents.size();
  CHECK(nprods == keys2prods.last());
  auto old_globals = old_mesh->ask_globals(ent_dim);
  auto comm = old_mesh->comm();
  auto old_ents2lins = copies_to_linear_owners(comm, old_globals);
  auto lins2old_ents = old_ents2lins.invert();
  auto nlins = lins2old_ents.nroots();
  auto lin_rep_counts = old_ents2lins.exch_reduce(global_rep_counts, 1, SUM);
  CHECK(lin_rep_counts.size() == nlins);
  auto lin_local_offsets = offset_scan(lin_rep_counts);
  auto lin_global_count = lin_local_offsets.last();
  auto lin_global_offset = comm->exscan(lin_global_count, SUM);
  Write<GO> lin_globals(nlins);
  auto write_lin_globals = LAMBDA(LO lin) {
    lin_globals[lin] = lin_local_offsets[lin] + lin_global_offset;
  };
  parallel_for(nlins, write_lin_globals);
  auto old_ents2new_globals = lins2old_ents.exch(
      Read<GO>(lin_globals), 1);
  Read<GO> same_ents2new_globals;
  Read<GO> prods2new_globals;
  auto edge2rep_order = LOs();
  if (ent_dim == VERT && key_dim == EDGE) {
    edge2rep_order = old_mesh->get_array<LO>(EDGE, "edge2rep_order");
  }
  find_new_offsets(old_ents2new_globals, same_ents2old_ents,
      keys2kds, keys2reps, keys2prods, edge2rep_order,
      &same_ents2new_globals, &prods2new_globals);
  auto nnew_ents = new_mesh->nents(ent_dim);
  CHECK(nnew_ents == nsame_ents + nprods);
  Write<GO> new_globals(nnew_ents);
  map_into(same_ents2new_globals, same_ents2new_ents, new_globals, 1);
  map_into(prods2new_globals, prods2new_ents, new_globals, 1);
  new_mesh->add_tag(ent_dim, "global", 1, OSH_GLOBAL, Read<GO>(new_globals));
}

static void modify_globals_cheap(Mesh* old_mesh, Mesh* new_mesh,
    Int ent_dim,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents) {
  auto nnew_ents = new_mesh->nents(ent_dim);
  if (!new_mesh->could_be_shared(ent_dim)) {
    auto start = new_mesh->comm()->exscan(GO(nnew_ents), SUM);
    auto globals = Read<GO>(nnew_ents, start, 1);
    new_mesh->add_tag(ent_dim, "global", 1, OSH_GLOBAL, globals);
  }
  auto old_owned = old_mesh->owned(ent_dim);
  auto same_owned = unmap(same_ents2old_ents, old_owned, 1);
  auto new_owned_w = Write<I8>(new_mesh->nents(ent_dim), 0);
  map_into(same_owned, same_ents2new_ents, new_owned_w, 1);
  auto nprods = prods2new_ents.size();
  map_into(Read<I8>(nprods, 1), prods2new_ents, new_owned_w, 1);
  auto new_owned = Read<I8>(new_owned_w);
  auto local_offsets = offset_scan(new_owned);
  auto nnew_owned = local_offsets.last();
  auto start = new_mesh->comm()->exscan(GO(nnew_owned), SUM);
  auto new_globals_w = Write<GO>(nnew_ents);
  parallel_for(nnew_ents, LAMBDA(LO e) {
    new_globals_w[e] = local_offsets[e] + start;
  });
  auto new_globals = Read<GO>(new_globals_w);
  new_globals = new_mesh->sync_array(ent_dim, new_globals, 1);
  new_mesh->add_tag(ent_dim, "global", 1, OSH_GLOBAL, new_globals);
}

void modify_ents(Mesh* old_mesh, Mesh* new_mesh,
    Int ent_dim, Int key_dim,
    LOs keys2kds,
    LOs keys2prods,
    LOs prod_verts2verts,
    LOs old_lows2new_lows,
    LOs* p_prods2new_ents,
    LOs* p_same_ents2old_ents,
    LOs* p_same_ents2new_ents,
    LOs* p_old_ents2new_ents) {
  *p_same_ents2old_ents = collect_same(old_mesh, ent_dim, key_dim,
      keys2kds);
  auto nkeys = keys2kds.size();
  CHECK(nkeys == keys2prods.size() - 1);
  auto keys2nprods = get_degrees(keys2prods);
  auto keys2reps = get_keys2reps(old_mesh, ent_dim, key_dim,
      keys2kds, keys2nprods);
  auto local_rep_counts = get_rep_counts(old_mesh, ent_dim,
      keys2reps, keys2nprods, *p_same_ents2old_ents, false);
  auto local_offsets = offset_scan(local_rep_counts);
  auto nnew_ents = local_offsets.last();
  auto edge2rep_order = LOs();
  if (ent_dim == VERT && key_dim == EDGE) {
    /* recompute this because the local version differs
       from the global one */
    auto edges_are_keys = mark_image(keys2kds, old_mesh->nents(EDGE));
    edge2rep_order = get_edge2rep_order(old_mesh, edges_are_keys);
  }
  find_new_offsets(local_offsets, *p_same_ents2old_ents,
      keys2kds, keys2reps, keys2prods, edge2rep_order,
      p_same_ents2new_ents, p_prods2new_ents);
  auto nold_ents = old_mesh->nents(ent_dim);
  *p_old_ents2new_ents = map_onto(*p_same_ents2new_ents, *p_same_ents2old_ents,
      nold_ents, -1, 1);
  if (ent_dim == VERT) {
    new_mesh->set_verts(nnew_ents);
  } else {
    modify_conn(old_mesh, new_mesh, ent_dim, prod_verts2verts,
        *p_prods2new_ents, *p_same_ents2old_ents, *p_same_ents2new_ents,
        old_lows2new_lows);
  }
  if (old_mesh->comm()->size() > 1) {
    modify_owners(old_mesh, new_mesh, ent_dim, *p_prods2new_ents,
        *p_same_ents2old_ents, *p_same_ents2new_ents, *p_old_ents2new_ents);
  }
  if (old_mesh->keeps_canonical_globals()) {
    auto global_rep_counts = get_rep_counts(old_mesh, ent_dim,
        keys2reps, keys2nprods, *p_same_ents2old_ents, true);
    modify_globals(old_mesh, new_mesh, ent_dim, key_dim,
        keys2kds, keys2prods, *p_prods2new_ents,
        *p_same_ents2old_ents, *p_same_ents2new_ents,
        keys2reps, global_rep_counts);
  } else {
    modify_globals_cheap(old_mesh, new_mesh, ent_dim,
        *p_prods2new_ents, *p_same_ents2old_ents, *p_same_ents2new_ents);
  }
}

void set_owners_by_indset(Mesh* mesh, Int key_dim, LOs keys2kds) {
  if (mesh->comm()->size() == 1) return;
  auto kd_owners = mesh->ask_owners(key_dim);
  auto nkeys = keys2kds.size();
  auto elem_dim = mesh->dim();
  auto kds2elems = mesh->ask_up(key_dim, elem_dim);
  auto kds2kd_elems = kds2elems.a2ab;
  auto kd_elems2elems = kds2elems.ab2b;
  auto elem_owners = mesh->ask_owners(elem_dim);
  auto elems2owners = mesh->ask_dist(elem_dim);
  auto new_own_ranks = deep_copy(elem_owners.ranks);
  auto f = LAMBDA(LO key) {
    auto kd = keys2kds[key];
    auto kd_rank = kd_owners.ranks[kd];
    for (auto kd_elem = kds2kd_elems[kd];
         kd_elem < kds2kd_elems[kd + 1];
         ++kd_elem) {
      auto elem = kd_elems2elems[kd_elem];
      new_own_ranks[elem] = kd_rank;
    }
  };
  parallel_for(nkeys, f);
  elem_owners = update_ownership(elems2owners, new_own_ranks);
  mesh->set_owners(elem_dim, elem_owners);
}
