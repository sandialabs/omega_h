void classify_sides_by_exposure(Mesh& mesh, Read<I8> side_is_exposed) {
  auto dim = mesh.dim();
  auto ns = mesh.nents(dim - 1);
  Write<I8> class_dim(ns);
  auto f = LAMBDA(LO s) {
    class_dim[s] = static_cast<I8>(dim - side_is_exposed[s]);
  };
  parallel_for(ns, f);
  mesh.add_tag<I8>(dim - 1, "class_dim", 1, class_dim);
}

void classify_hinges_by_sharpness(Mesh& mesh,
    Read<I8> hinge_is_exposed,
    Read<I8> hinge_is_sharp) {
  auto dim = mesh.dim();
  auto nh = mesh.nents(dim - 2);
  Write<I8> class_dim(nh);
  auto f = LAMBDA(LO h) {
    class_dim[h] = static_cast<I8>(dim - hinge_is_exposed[h] - hinge_is_sharp[h]);
  };
  parallel_for(nh, f);
  mesh.add_tag<I8>(dim - 2, "class_dim", 1, class_dim);
}

void classify_vertices_by_sharp_edges(Mesh& mesh,
    Read<I8> vert_is_exposed,
    Read<I8> edge_is_sharp) {
  auto nv = mesh.nents(VERT);
  auto v2e = mesh.ask_adj(VERT,EDGE);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  Write<I8> class_dim(nv);
  auto f = LAMBDA(LO v) {
    auto begin = v2ve[v];
    auto end = v2ve[v + 1];
    Int nadj_sharp = 0;
    for (auto ve = begin; ve < end; ++ve) {
      auto e = ve2e[ve];
      nadj_sharp += edge_is_sharp[e];
    }
    if (nadj_sharp == 0) {
      class_dim[v] = 3 - vert_is_exposed[v];
    } else if (nadj_sharp == 2) {
      class_dim[v] = 1;
    } else {
      class_dim[v] = 0;
    }
  };
  parallel_for(nv, f);
  mesh.add_tag<I8>(VERT, "class_dim", 1, class_dim);
}

void classify_by_angles(Mesh& mesh, Real sharp_angle) {
  auto dim = mesh.dim();
  mesh.add_tag<I8>(dim, "class_dim", 1, Read<I8>(mesh.nelems(), static_cast<I8>(dim)));
  auto side_is_exposed = mark_exposed_sides(mesh);
  classify_sides_by_exposure(mesh, side_is_exposed);
  auto hinge_is_exposed = mark_down(mesh, dim - 1, dim - 2, side_is_exposed);
  auto surf_side2side = collect_marked(side_is_exposed);
  auto surf_side_normals = surf::get_side_normals(mesh, surf_side2side);
  auto surf_hinge2hinge = collect_marked(hinge_is_exposed);
  auto nsurf_hinges = surf_hinge2hinge.size();
  auto nsides = mesh.nents(dim - 1);
  auto side2surf_side = invert_injective_map(surf_side2side, nsides);
  auto surf_hinge_angles = surf::get_hinge_angles(mesh,
      surf_side_normals, surf_hinge2hinge, side2surf_side);
  auto nhinges = mesh.nents(dim - 2);
  Write<I8> hinge_is_sharp(nhinges, 0);
  auto f = LAMBDA(LO surf_hinge) {
    LO hinge = surf_hinge2hinge[surf_hinge];
    hinge_is_sharp[hinge] = (surf_hinge_angles[surf_hinge] >= sharp_angle);
  };
  parallel_for(nsurf_hinges, f);
  classify_hinges_by_sharpness(mesh, hinge_is_exposed, hinge_is_sharp);
  if (dim == 2)
    return;
  auto vert_is_exposed = mark_down(mesh, 2, 0, side_is_exposed);
  classify_vertices_by_sharp_edges(mesh,
      vert_is_exposed, hinge_is_sharp);
}
