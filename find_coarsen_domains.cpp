static void mark_dead_ents(Mesh& mesh,
    LOs rails2edges,
    Read<I8> rail_col_verts,
    Int cell_dim,
    Write<I8>& dead_cells,
    Write<I8>& dead_sides) {
  auto e2c = mesh.ask_up(EDGE, cell_dim);
  auto e2ec = e2c.a2ab;
  auto ec2c = e2c.ab2b;
  auto ec_codes = e2c.codes;
  auto cs2s = mesh.ask_down(cell_dim, cell_dim - 1).ab2b;
  auto nccs = simplex_degrees[cell_dim][cell_dim - 1];
  auto nrails = rails2edges.size();
  auto f = LAMBDA(LO rail) {
    auto e = rails2edges[rail];
    auto eev_col = rail_col_verts[rail];
    auto eev_onto = 1 - eev_col;
    for (auto ec = e2ec[e]; ec < e2ec[e + 1]; ++ec) {
      auto c = ec2c[ec];
      auto ec_code = ec_codes[ec];
      auto cce = code_which_down(ec_code);
      auto rot = code_rotation(ec_code);
      auto cev_onto = rot ^ eev_onto;
      auto ccv_onto = down_templates[cell_dim][EDGE][cce][cev_onto];
      auto ccs_opp = opposite_templates[cell_dim][VERT][ccv_onto];
      auto s_opp = cs2s[c * nccs + ccs_opp];
      dead_cells[c] = 1;
      dead_sides[s_opp] = 1;
    }
  };
  parallel_for(nrails, f);
}

Few<Read<I8>, 4> mark_dead_ents(Mesh& mesh,
    LOs rails2edges,
    Read<I8> rail_col_verts) {
  Few<Write<I8>, 4> writes;
  writes[EDGE] = deep_copy(mark_image(rails2edges, mesh.nedges()));
  for (Int dim = EDGE + 1; dim <= mesh.dim(); ++dim)
    writes[dim] = Write<I8>(mesh.nents(dim), 0);
  for (Int dim = mesh.dim(); dim > EDGE; --dim)
    mark_dead_ents(mesh, rails2edges, rail_col_verts,
        dim, writes[dim], writes[dim - 1]);
  Few<Read<I8>, 4> reads;
  for (Int dim = 0; dim < 4; ++dim)
    reads[dim] = writes[dim];
  return reads;
}

Adj find_coarsen_domains(Mesh& mesh,
    LOs keys2verts,
    Int ent_dim,
    Read<I8> ents_are_dead) {
  auto nkeys = keys2verts.size();
  auto v2e = mesh.ask_up(VERT, ent_dim);
  auto k2e = unmap_adjacency(keys2verts, v2e);
  auto k2ke = k2e.a2ab;
  auto ke2e = k2e.ab2b;
  auto ke_codes = k2e.codes;
  auto ke2k = invert_fan(k2ke);
  auto ents_are_live = invert_marks(ents_are_dead);
  auto kes_are_live = unmap(ke2e, ents_are_live, 1);
  auto lke2ke = collect_marked(kes_are_live);
  auto lke2k = unmap(lke2ke, ke2k, 1);
  auto lke_codes = unmap(lke2ke, ke_codes, 1);
  auto lke2e = unmap(lke2ke, ke2e, 1);
  auto k2lke = invert_funnel(lke2k, nkeys);
  return Adj(k2lke, lke2e, lke_codes);
}
