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
