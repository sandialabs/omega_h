/* Construct a graph from each vertex
   to the entities which use that vertex
   as their first vertex in the canonical
   downward ordering.
   The point of this is just to establish,
   for each entity, a single vertex that is
   "responsible for" that entity */
Graph find_entities_of_first_vertices(
    Mesh* mesh, Int ent_dim) {
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
  auto fv2fve = offset_scan(LOs(degrees));
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

LOs ent_order_from_vert_order(Mesh* mesh,
    Int ent_dim, LOs new_verts2old_verts) {
  auto old_verts2old_ents = find_entities_of_first_vertices(mesh, ent_dim);
  auto new_verts2old_ents = unmap_graph(new_verts2old_verts, old_verts2old_ents);
  auto new_ents2old_ents = new_verts2old_ents.ab2b;
  CHECK(new_ents2old_ents.size() == mesh.nents(ent_dim));
  return new_ents2old_ents;
}

void reorder_mesh(Mesh* old_mesh, Mesh* new_mesh,
    LOs new_verts2old_verts) {
  LOs new_ents2old_ents[4];
  new_ents2old_ents[VERT] = new_verts2old_verts;
  for (Int ent_dim = 1; ent_dim <= old_mesh.dim(); ++ent_dim) {
    new_ents2old_ents[ent_dim] = ent_order_from_vert_order(old_mesh,
        ent_dim, new_verts2old_verts);
  }
  unmap_mesh(old_mesh, new_mesh, new_ents2old_ents);
}

void reorder_mesh(Mesh* mesh,
    LOs new_verts2old_verts) {
  Mesh new_mesh;
  reorder_mesh(mesh, new_mesh, new_verts2old_verts);
  mesh = new_mesh;
}

void reorder_by_hilbert(Mesh* mesh) {
  auto coords = mesh.coords();
  LOs new_verts2old_verts = hilbert::sort_coords(coords, mesh.dim());
  reorder_mesh(mesh, new_verts2old_verts);
}
