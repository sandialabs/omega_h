static void add_ents2verts(Mesh& mesh, Int edim, LOs ev2v) {
  if (edim == 1) {
    mesh.set_ents(edim, Adj(ev2v));
  } else {
    Int ldim = edim - 1;
    LOs lv2v = mesh.ask_adj(ldim, VERT).ab2b;
    Adj v2l = mesh.ask_adj(VERT, ldim);
    Adj down = reflect_down(ev2v, lv2v, v2l, edim, ldim);
    mesh.set_ents(edim, down);
  }
}

void build_from_elems2verts(Mesh& mesh, Int edim, LOs ev2v, LO nverts) {
  mesh.set_dim(edim);
  mesh.set_verts(nverts);
  for (Int mdim = 1; mdim < edim; ++mdim) {
    LOs mv2v = find_unique(ev2v, edim, mdim);
    add_ents2verts(mesh, mdim, mv2v);
  }
  add_ents2verts(mesh, edim, ev2v);
}
