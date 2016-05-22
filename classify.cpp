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
