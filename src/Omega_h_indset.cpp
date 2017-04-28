#include "Omega_h_indset.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_mesh.hpp"

namespace Omega_h {

namespace indset {

enum { NOT_IN, IN, UNKNOWN };

static Read<I8> local_iteration(
    LOs xadj, LOs adj, Reals quality, Read<GO> global, Read<I8> old_state) {
  auto n = global.size();
  Write<I8> new_state = deep_copy(old_state);
  auto f = OMEGA_H_LAMBDA(LO v) {
    if (old_state[v] != UNKNOWN) return;
    auto begin = xadj[v];
    auto end = xadj[v + 1];
    // nodes adjacent to chosen ones are rejected
    for (auto j = begin; j < end; ++j) {
      auto u = adj[j];
      if (old_state[u] == IN) {
        new_state[v] = NOT_IN;
        return;
      }
    }
    // check if node is a local maximum
    auto v_qual = quality[v];
    for (auto j = begin; j < end; ++j) {
      auto u = adj[j];
      // neighbor was rejected, ignore its presence
      if (old_state[u] == NOT_IN) continue;
      auto u_qual = quality[u];
      // neighbor has higher quality
      if (u_qual > v_qual) return;
      // neighbor has equal quality, tiebreaker by global ID
      if (u_qual == v_qual && global[u] > global[v]) return;
    }
    // only local maxima reach this line
    new_state[v] = IN;
  };
  parallel_for(n, f);
  return new_state;
}

static Read<I8> iteration(Mesh* mesh, Int dim, LOs xadj, LOs adj, Reals quality,
    Read<GO> global, Read<I8> old_state) {
  auto local_state = local_iteration(xadj, adj, quality, global, old_state);
  auto synced_state = mesh->sync_array(dim, local_state, 1);
  return synced_state;
}

static Read<I8> find(Mesh* mesh, Int dim, LOs xadj, LOs adj, Reals quality,
    Read<GO> global, Read<I8> candidates) {
  auto n = global.size();
  OMEGA_H_CHECK(quality.size() == n);
  OMEGA_H_CHECK(candidates.size() == n);
  auto initial_state = Write<I8>(n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    if (candidates[i])
      initial_state[i] = UNKNOWN;
    else
      initial_state[i] = NOT_IN;
  };
  parallel_for(n, f);
  auto comm = mesh->comm();
  auto state = Read<I8>(initial_state);
  while (get_max(comm, state) == UNKNOWN) {
    state = iteration(mesh, dim, xadj, adj, quality, global, state);
  }
  return state;
}
}  // namespace indset

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Graph graph, Reals quality, Read<I8> candidates) {
  auto xadj = graph.a2ab;
  auto adj = graph.ab2b;
  auto globals = mesh->ask_globals(ent_dim);
  return indset::find(mesh, ent_dim, xadj, adj, quality, globals, candidates);
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
