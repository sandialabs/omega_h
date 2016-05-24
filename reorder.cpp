/* Construct a graph from each vertex
   to the entities which use that vertex
   as their first vertex in the canonical
   downward ordering.
   The point of this is just to establish,
   for each entity, a single vertex that is
   "responsible for" that entity */
Graph find_entities_of_first_vertices(
    Mesh& mesh, Int ent_dim) {
  auto v2e = mesh.ask_up(VERT, ent_dim);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ve_codes = v2e.codes;
  auto nv = mesh.nverts();
  Write<LO> degrees(nv);
  auto count = LAMBDA(LO v) {
    Int n = 0;
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      n += (code_which_down(ve_codes[ve]) == 0);
    }
    degrees[v] = n;
  };
  parallel_for(nv, count);
  auto fv2fve = offset_scan<LO,LO>(degrees);
  Write<LO> fve2e(fv2fve.last());
  auto fill = LAMBDA(LO v) {
    LO fve = fv2fve[v];
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      if (code_which_down(ve_codes[ve]) == 0) {
        fve2e[fve++] = ve2e[ve];
      }
    }
  };
  parallel_for(nv, fill);
  return Graph(fv2fve, fve2e);
}

LOs ent_order_from_vert_order(Mesh& mesh,
    Int ent_dim, LOs new_vert2old_vert) {
  auto old_verts2old_ents = find_entities_of_first_vertices(mesh, ent_dim);
  auto new_verts2old_ents = unmap_graph(new_vert2old_vert, old_verts2old_ents);
  auto new_ents2old_ents = new_verts2old_ents.ab2b;
  CHECK(new_ents2old_ents.size() == mesh.nents(ent_dim));
  return new_ents2old_ents;
}

void reorder_tags(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim, LOs new_ents2old_ents) {
  for (Int i = 0; i < old_mesh.ntags(ent_dim); ++i) {
    auto tag = old_mesh.get_tag(ent_dim, i);
    if (is<I8>(tag)) {
      new_mesh.add_tag<I8>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, to<I8>(tag)->array(), tag->ncomps()));
    } else if (is<I32>(tag)) {
      new_mesh.add_tag<I32>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, to<I32>(tag)->array(), tag->ncomps()));
    } else if (is<I64>(tag)) {
      new_mesh.add_tag<I64>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, to<I64>(tag)->array(), tag->ncomps()));
    } else if (is<Real>(tag)) {
      new_mesh.add_tag<Real>(ent_dim, tag->name(), tag->ncomps(),
          unmap(new_ents2old_ents, to<Real>(tag)->array(), tag->ncomps()));
    }
  }
}

void reorder_down(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim, LOs new_ents2old_ents,
    LOs old_lows2new_lows) {
  auto deg = simplex_degrees[ent_dim][ent_dim - 1];
  auto old_ents2old_lows = old_mesh.ask_down(ent_dim, ent_dim - 1);
  auto oel2ol = old_ents2old_lows.ab2b;
  auto oe2l_codes = old_ents2old_lows.codes;
  auto nel2ol = unmap(new_ents2old_ents, oel2ol, deg);
  auto nel2nl = compound_maps(nel2ol, old_lows2new_lows);
  auto new_ents2new_lows = Adj(nel2nl);
  if (oe2l_codes.size()) {
    auto ne2l_codes = unmap(new_ents2old_ents, oe2l_codes, deg);
    new_ents2new_lows.codes = ne2l_codes;
  }
  new_mesh.set_ents(ent_dim, new_ents2new_lows);
}

void reorder_mesh(Mesh& old_mesh, Mesh& new_mesh,
    LOs new_vert2old_vert) {
  new_mesh.set_dim(old_mesh.dim());
  new_mesh.set_verts(old_mesh.nverts());
  reorder_tags(old_mesh, new_mesh, VERT, new_vert2old_vert);
  auto old_lows2new_lows = invert_permutation(new_vert2old_vert);
  for (Int ent_dim = 1; ent_dim <= old_mesh.dim(); ++ent_dim) {
    auto new_ents2old_ents = ent_order_from_vert_order(old_mesh,
        ent_dim, new_vert2old_vert);
    reorder_down(old_mesh, new_mesh, ent_dim, new_ents2old_ents,
        old_lows2new_lows);
    reorder_tags(old_mesh, new_mesh, ent_dim, new_ents2old_ents);
    old_lows2new_lows = invert_permutation(new_ents2old_ents);
  }
}

void reorder_mesh(Mesh& mesh,
    LOs new_verts2old_verts) {
  Mesh new_mesh;
  reorder_mesh(mesh, new_mesh, new_verts2old_verts);
  mesh = new_mesh;
}
