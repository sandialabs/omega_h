#include "Omega_h_swap.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_profile.hpp"
#include "Omega_h_swap2d.hpp"
#include "Omega_h_swap3d.hpp"

#include <iostream>

namespace Omega_h {

bool swap_part1(Mesh* mesh, AdaptOpts const& opts) {
  mesh->set_parting(OMEGA_H_GHOSTED);
  auto comm = mesh->comm();
  auto elems_are_cands =
      mark_sliver_layers(mesh, opts.min_quality_desired, opts.nsliver_layers);
  OMEGA_H_CHECK(get_max(comm, elems_are_cands) == 1);
  auto edges_are_cands = mark_down(mesh, mesh->dim(), EDGE, elems_are_cands);
  /* only swap interior edges */
  auto edges_are_inter = mark_by_class_dim(mesh, EDGE, mesh->dim());
  edges_are_cands = land_each(edges_are_cands, edges_are_inter);
  if (get_max(comm, edges_are_cands) <= 0) return false;
  mesh->add_tag(EDGE, "candidate", 1, edges_are_cands);
  return true;
}

void filter_swap(Read<I8> keep_cands, LOs* cands2edges, Reals* cand_quals) {
  auto kept2old = collect_marked(keep_cands);
  *cands2edges = unmap(kept2old, *cands2edges, 1);
  if (cand_quals) *cand_quals = unmap(kept2old, *cand_quals, 1);
}

Read<I8> filter_swap_improve(Mesh* mesh, LOs cands2edges, Reals cand_quals) {
  OMEGA_H_CHECK(mesh->owners_have_all_upward(EDGE));
  auto elem_quals = mesh->ask_qualities();
  auto edges2elems = mesh->ask_up(EDGE, mesh->dim());
  auto edge_old_quals = graph_reduce(edges2elems, elem_quals, 1, OMEGA_H_MIN);
  edge_old_quals = mesh->sync_array(EDGE, edge_old_quals, 1);
  auto cand_old_quals = read(unmap(cands2edges, edge_old_quals, 1));
  return gt_each(cand_quals, cand_old_quals);
}

bool swap_edges(Mesh* mesh, AdaptOpts const& opts) {
  OMEGA_H_TIME_FUNCTION;
  bool ret = false;
  if (mesh->dim() == 3)
    ret = swap_edges_3d(mesh, opts);
  else if (mesh->dim() == 2)
    ret = swap_edges_2d(mesh, opts);
  return ret;
}

}  // end namespace Omega_h
