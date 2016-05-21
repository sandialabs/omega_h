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

void build_from_elems_and_coords(Mesh& mesh, Int edim, LOs ev2v, Reals coords) {
  LO nverts = coords.size() / edim;
  build_from_elems2verts(mesh, edim, ev2v, nverts);
  mesh.add_coords();
  mesh.set_coords(coords);
}

void build_box(Mesh& mesh,
    Real x, Real y, Real z,
    LO nx, LO ny, LO nz) {
  CHECK(nx > 0);
  CHECK(ny > 0);
  CHECK(nz >= 0);
  if (nz == 0) {
    LOs qv2v;
    Reals coords;
    make_2d_box(x, y, nx, ny, qv2v, coords);
    LOs tv2v = simplify::tris_from_quads(qv2v);
    build_from_elems_and_coords(mesh, TRI, tv2v, coords);
  } else {
    LOs hv2v;
    Reals coords;
    make_3d_box(x, y, z, nx, ny, nz, hv2v, coords);
    LOs tv2v = simplify::tets_from_hexes(hv2v);
    build_from_elems_and_coords(mesh, TET, tv2v, coords);
  }
}
