static bool refine_ghosted(Mesh& mesh, Real min_qual) {
  auto comm = mesh.comm();
  auto edges_are_cands = mesh.get_array<I8>(EDGE, "candidate");
  mesh.remove_tag(EDGE, "candidate");
  auto cands2edges = collect_marked(edges_are_cands);
  auto cand_quals = refine_qualities(mesh, cands2edges);
  auto cands_are_good = each_geq_to(cand_quals, min_qual);
  if (comm->allreduce(max(cands_are_good), MAX) == 0) return false;
  auto nedges = mesh.nents(EDGE);
  auto edges_are_initial_w = Write<I8>(nedges, 0);
  map_into(cands_are_good, cands2edges, edges_are_initial_w, 1);
  auto edges_are_initial = Read<I8>(edges_are_initial_w);
  auto edge_quals_w = Write<Real>(nedges, 0.);
  map_into(cand_quals, cands2edges, edge_quals_w, 1);
  auto edge_quals = Reals(edge_quals_w);
  auto edges_are_keys = find_indset(mesh, EDGE, edge_quals, edges_are_initial);
  mesh.add_tag(EDGE, "key", 1, OSH_DONT_TRANSFER, edges_are_keys);
  mesh.add_tag(EDGE, "edge2rep_order", 1, OSH_DONT_TRANSFER,
      get_edge2rep_order(mesh, edges_are_keys));
  auto keys2edges = collect_marked(edges_are_keys);
  set_owners_by_indset(mesh, EDGE, keys2edges);
  return true;
}

static void refine_element_based(Mesh& mesh) {
  auto comm = mesh.comm();
  auto edges_are_keys = mesh.get_array<I8>(EDGE, "key");
  auto keys2edges = collect_marked(edges_are_keys);
  auto nkeys = keys2edges.size();
  auto new_mesh = Mesh();
  new_mesh.set_comm(comm);
  new_mesh.set_dim(mesh.dim());
  new_mesh.set_partition(mesh.partition());
  auto keys2midverts = LOs();
  auto old_verts2new_verts = LOs();
  auto old_lows2new_lows = LOs();
  for (Int ent_dim = 0; ent_dim <= mesh.dim(); ++ent_dim) {
    auto keys2prods = LOs();
    auto prod_verts2verts = LOs();
    if (ent_dim == VERT) {
      keys2prods = LOs(nkeys + 1, 0, 1);
    } else {
      refine_products(mesh, ent_dim, keys2edges, keys2midverts,
          old_verts2new_verts, keys2prods, prod_verts2verts);
    }
    auto prods2new_ents = LOs();
    auto same_ents2old_ents = LOs();
    auto same_ents2new_ents = LOs();
    auto old_ents2new_ents = LOs();
    modify_ents(mesh, new_mesh, ent_dim, EDGE, keys2edges,
        keys2prods, prod_verts2verts, old_lows2new_lows,
        prods2new_ents,
        same_ents2old_ents, same_ents2new_ents,
        old_ents2new_ents);
    if (ent_dim == VERT) {
      keys2midverts = prods2new_ents;
      old_verts2new_verts = old_ents2new_ents;
      transfer_linear_interp(mesh, new_mesh, keys2edges, keys2midverts,
          same_ents2new_ents);
      transfer_metric(mesh, new_mesh, keys2edges, keys2midverts,
          same_ents2new_ents);
    }
    old_lows2new_lows = old_ents2new_ents;
  }
  mesh = new_mesh;
}

bool refine(Mesh& mesh, Real min_qual) {
  mesh.set_partition(GHOSTED);
  if (!refine_ghosted(mesh, min_qual)) return false;
  mesh.set_partition(ELEMENT_BASED);
  refine_element_based(mesh);
  return true;
}
