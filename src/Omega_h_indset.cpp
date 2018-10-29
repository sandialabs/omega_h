#include "Omega_h_indset_inline.hpp"

namespace Omega_h {

struct QualityCompare {
  Reals quality;
  GOs global;
  OMEGA_H_DEVICE bool operator()(LO u, LO v) const {
    auto const v_qual = quality[v];
    auto const u_qual = quality[u];
    if (u_qual != v_qual) return u_qual < v_qual;
    // neighbor has equal quality, tiebreaker by global ID
    return global[u] < global[v];
  }
};

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Graph graph, Reals quality, Read<I8> candidates) {
  auto xadj = graph.a2ab;
  auto adj = graph.ab2b;
  QualityCompare compare;
  compare.quality = quality;
  compare.global = mesh->globals(ent_dim);
  return indset::find(mesh, ent_dim, xadj, adj, candidates, compare);
}

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Reals quality, Read<I8> candidates) {
  if (ent_dim == mesh->dim()) return candidates;
  mesh->owners_have_all_upward(ent_dim);
  OMEGA_H_CHECK(mesh->owners_have_all_upward(ent_dim));
  auto graph = mesh->ask_star(ent_dim);
  return find_indset(mesh, ent_dim, graph, quality, candidates);
}

}  // end namespace Omega_h
