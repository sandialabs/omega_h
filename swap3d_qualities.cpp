template <typename Measure>
static void swap3d_qualities_tmpl(Mesh* mesh, LOs cands2edges,
    Reals* cand_quals, Read<I8>* cand_configs) {
  auto edges2tets = mesh->ask_up(EDGE, TET);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edge_tets2tets = edges2tets.ab2b;
  auto edge_tet_codes = edges2tets.codes;
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto tet_verts2verts = mesh->ask_verts_of(TET);
  auto edges_are_owned = mesh->owned(EDGE);
  Measure measure(mesh);
  auto ncands = cands2edges.size();
  auto cand_quals_w = Write<Real>(ncands);
  auto cand_configs_w = Write<I8>(ncands);
  auto f = LAMBDA(LO cand) {
    auto edge = cands2edges[cand];
    /* non-owned edges will have incomplete cavities
       and will run into the topological assertions
       in find_loop(). don't bother; their results
       will be overwritten by the owner's anyways */
    if (!edges_are_owned[edge]) {
      cand_configs_w[cand] = -1;
      cand_quals_w[cand] = -1.0;
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
      cand_configs_w[cand] = -1;
      cand_quals_w[cand] = -1.0;
      return;
    }
    auto choice = swap3d::choose(loop, measure);
    static_assert(swap3d::MAX_CONFIGS <= INT8_MAX,
        "int8_t must be able to represent all swap configurations");
    cand_configs_w[cand] = static_cast<I8>(choice.mesh);
    cand_quals_w[cand] = choice.quality;
  };
  parallel_for(ncands, f);
  *cand_quals = cand_quals_w;
  *cand_configs = cand_configs_w;
}

void swap3d_qualities(Mesh* mesh, LOs cands2edges,
    Reals* cand_quals, Read<I8>* cand_configs) {
  CHECK(mesh->parting() == OSH_GHOSTED);
  CHECK(mesh->dim() == 3);
  if (mesh->has_tag(VERT, "metric")) {
    swap3d_qualities_tmpl<MetricElementQualities>(mesh, cands2edges,
        cand_quals, cand_configs);
  } else {
    swap3d_qualities_tmpl<RealElementQualities>(mesh, cands2edges,
        cand_quals, cand_configs);
  }
  *cand_quals = mesh->sync_subset_array(EDGE, *cand_quals, cands2edges, -1.0, 1);
  *cand_configs = mesh->sync_subset_array(EDGE, *cand_configs, cands2edges, I8(-1), 1);
}
