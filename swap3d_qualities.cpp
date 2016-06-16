template <typename Measure>
static Reals swap3d_qualities_tmpl(Mesh* mesh, LOs cands2edges) {
  auto edges2tets = mesh->ask_up(EDGE, TET);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edge_tets2tets = edges2tets.ab2b;
  auto edge_tet_codes = edges2tets.codes;
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto tet_verts2verts = mesh->ask_verts_of(TET);
  auto edges_are_owned = mesh->owned(EDGE);
  Measure measure(mesh);
  auto ncands = cands2edges.size();
  auto cand_quals = Write<Real>(ncands);
  auto f = LAMBDA(LO cand) {
    auto edge = cands2edges[cand];
    /* non-owned edges will have incomplete cavities
       and will run into the topological assertions
       in find_loop(). don't bother; their results
       will be overwritten by the owner's anyways */
    if (!edges_are_owned[edge]) {
      cand_quals[cand] = -1.0;
      return;
    }
    auto loop = swap3d::find_loop(
        edges2edge_tets,
        edge_tets2tets,
        edge_tet_codes,
        edge_verts2verts,
        tet_verts2verts,
        edge);
    if (loop.size > swap3d::MAX_EDGE_SWAP) {
      cand_quals[cand] = -1.0;
      return;
    }
    auto choice = swap3d::choose(loop, measure);
    if (choice.mesh == -1) {
      cand_quals[cand] = -1.0;
      return;
    }
    cand_quals[cand] = choice.quality;
  };
  parallel_for(ncands, f);
  return cand_quals;
}

Reals swap3d_qualities(Mesh* mesh, LOs cands2edges) {
  CHECK(mesh->partition() == GHOSTED);
  auto cand_quals = Reals();
  CHECK(mesh->dim() == 3);
  if (mesh->has_tag(VERT, "metric")) {
    cand_quals = swap3d_qualities_tmpl<MetricElementQualities>(mesh, cands2edges);
  } else {
    cand_quals = swap3d_qualities_tmpl<RealElementQualities>(mesh, cands2edges);
  }
  return mesh->sync_subset_array(EDGE, cand_quals, cands2edges, -1.0, 1);
}
