bool swap_part1(Mesh* mesh, Real qual_ceil, Int nlayers) {
  mesh->set_partition(GHOSTED);
  auto comm = mesh->comm();
  auto elems_are_cands = mark_sliver_layers(mesh, qual_ceil, nlayers);
  CHECK(comm->allreduce(max(elems_are_cands), OSH_MAX) == 1);
  auto edges_are_cands = mark_down(mesh, mesh->dim(), EDGE, elems_are_cands);
  auto edges_are_inter = mark_by_class_dim(mesh, EDGE, mesh->dim());
  edges_are_cands = land_each(edges_are_cands, edges_are_inter);
  /* only swap interior edges */
  if (comm->reduce_and(max(edges_are_cands) <= 0)) return false;
  mesh->add_tag(EDGE, "candidate", 1, OSH_DONT_TRANSFER, edges_are_cands);
  return true;
}

void filter_swap_improve(Mesh* mesh,
    LOs* cands2edges, Reals* cand_quals) {
  CHECK(mesh->owners_have_all_upward(EDGE));
  auto elem_quals = mesh->ask_qualities();
  auto edges2elems = mesh->ask_up(EDGE, mesh->dim());
  auto edge_old_quals = graph_reduce(edges2elems, elem_quals, 1, OSH_MIN);
  edge_old_quals = mesh->sync_array(EDGE, edge_old_quals, 1);
  auto cand_old_quals = unmap(*cands2edges, edge_old_quals, 1);
  auto keep_cands = gt_each(*cand_quals, cand_old_quals);
  auto kept2old = collect_marked(keep_cands);
  *cands2edges = unmap(kept2old, *cands2edges, 1);
  *cand_quals = unmap(kept2old, *cand_quals, 1);
}

bool swap_edges(Mesh* mesh, Real qual_ceil, Int nlayers) {
  if (mesh->dim() == 3) return run_swap3d(mesh, qual_ceil, nlayers);
  if (mesh->dim() == 2) return swap2d(mesh, qual_ceil, nlayers);
  return false;
}
