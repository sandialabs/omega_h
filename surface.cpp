namespace surf {

namespace {

Reals get_triangle_normals(Mesh* mesh, LOs surf_tri2tri) {
  CHECK(mesh->dim() == 3);
  auto nsurf_tris = surf_tri2tri.size();
  auto fv2v = mesh->ask_verts_of(TRI);
  auto coords = mesh->coords();
  Write<Real> normals(nsurf_tris * 3);
  auto lambda = LAMBDA(LO surf_tri) {
    auto f = surf_tri2tri[surf_tri];
    auto v = gather_verts<3>(fv2v, f);
    auto x = gather_vectors<3, 3>(coords, v);
    auto b = simplex_basis<3, 2>(x);
    auto n = normalize(cross(b[0], b[1]));
    set_vec(normals, surf_tri, n);
  };
  parallel_for(nsurf_tris, lambda);
  return normals;
}

Reals get_edge_normals(Mesh* mesh, LOs surf_edge2edge) {
  CHECK(mesh->dim() == 2);
  auto nsurf_edges = surf_edge2edge.size();
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto coords = mesh->coords();
  Write<Real> normals(nsurf_edges * 2);
  auto lambda = LAMBDA(LO surf_edge) {
    auto e = surf_edge2edge[surf_edge];
    auto v = gather_verts<2>(ev2v, e);
    auto x = gather_vectors<2, 2>(coords, v);
    auto b = simplex_basis<2, 1>(x);
    auto n = normalize(perp(b[0]));
    set_vec(normals, surf_edge, n);
  };
  parallel_for(nsurf_edges, lambda);
  return normals;
}

template <Int dim>
Reals get_hinge_angles_tmpl(Mesh* mesh,
    Reals surf_side_normals,
    LOs surf_hinge2hinge,
    LOs side2surf_side) {
  auto nsurf_hinges = surf_hinge2hinge.size();
  auto hinges2sides = mesh->ask_up(dim - 2, dim - 1);
  auto hinges2hinge_sides = hinges2sides.a2ab;
  auto hinge_sides2sides = hinges2sides.ab2b;
  Write<Real> angles(nsurf_hinges);
  auto f = LAMBDA(LO surf_hinge) {
    auto hinge = surf_hinge2hinge[surf_hinge];
    auto begin = hinges2hinge_sides[hinge];
    auto end = hinges2hinge_sides[hinge + 1];
    Int i = 0;
    Vector<dim> n[2];
    for (auto hs = begin; hs < end; ++hs) {
      auto s = hinge_sides2sides[hs];
      auto ss = side2surf_side[s];
      if (-1 == ss) continue;
      n[i++] = get_vec<dim>(surf_side_normals, ss);
    }
    angles[surf_hinge] = acos(n[0] * n[1]);
  };
  parallel_for(nsurf_hinges, f);
  return angles;
}

} //end anonymous namespace

Reals get_side_normals(Mesh* mesh, LOs surf_side2side) {
  if (mesh->dim() == 3) {
    return get_triangle_normals(mesh, surf_side2side);
  } else {
    return get_edge_normals(mesh, surf_side2side);
  }
}

Reals get_hinge_angles(Mesh* mesh,
    Reals surf_side_normals,
    LOs surf_hinge2hinge,
    LOs side2surf_side) {
  if (mesh->dim() == 3) {
    return get_hinge_angles_tmpl<3>(mesh,
        surf_side_normals, surf_hinge2hinge, side2surf_side);
  } else {
    return get_hinge_angles_tmpl<2>(mesh,
        surf_side_normals, surf_hinge2hinge, side2surf_side);
  }
}

} //end namespace surface
