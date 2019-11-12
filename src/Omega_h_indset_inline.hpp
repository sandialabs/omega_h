#ifndef OMEGA_H_INDSET_INLINE_HPP
#define OMEGA_H_INDSET_INLINE_HPP

#include <Omega_h_array_ops.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_indset.hpp>
#include <Omega_h_mesh.hpp>

namespace Omega_h {
namespace indset {

enum { NOT_IN, IN, UNKNOWN };

template <class Compare>
inline Read<I8> local_iteration(
    LOs xadj, LOs adj, Read<I8> old_state, Compare compare) {
  auto n = xadj.size() - 1;
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
    for (auto j = begin; j < end; ++j) {
      auto u = adj[j];
      // neighbor was rejected, ignore its presence
      if (old_state[u] == NOT_IN) continue;
      if (!compare(u, v)) return;
    }
    // only local maxima reach this line
    new_state[v] = IN;
  };
  parallel_for(n, std::move(f));
  return new_state;
}

template <class Compare>
Read<I8> iteration(Mesh* mesh, Int dim, LOs xadj, LOs adj, Read<I8> old_state,
    Compare compare) {
  auto local_state = local_iteration(xadj, adj, old_state, compare);
  auto synced_state = mesh->sync_array(dim, local_state, 1);
  return synced_state;
}

template <class Compare>
Read<I8> find(Mesh* mesh, Int dim, LOs xadj, LOs adj, Read<I8> candidates,
    Compare compare) {
  auto n = xadj.size() - 1;
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
    state = iteration(mesh, dim, xadj, adj, state, compare);
  }
  return state;
}
}  // namespace indset
}  // namespace Omega_h

#endif
