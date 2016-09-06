#include "internal.hpp"
#include "loop.hpp"
#include "scan.hpp"
#include "simplices.hpp"
#include "indset.hpp"

namespace Omega_h {

namespace {

DEVICE Int add_unique(LO* stack, Int n, LO e, Int stack_max) {
  for (Int i = 0; i < n; ++i)
    if (stack[i] == e) return n;
  CHECK(n < stack_max);
  stack[n++] = e;
  return n;
}

DEVICE Int find_distance_2_elems(
    Graph const& keys2elems,
    LOs const& elems2verts,
    Int nverts_per_elem,
    Graph const& verts2elems,
    Read<I8> const& indset,
    LO key,
    LO* stack,
    Int stack_max) {
  if (!indset[key]) return 0;
  Int n = 0;
  for (auto key_elem = keys2elems.a2ab[key];
       key_elem < keys2elems.a2ab[key + 1];
       ++key_elem) {
    auto elem1 = keys2elems.ab2b[key_elem];
    for (auto elem_vert = elem1 * nverts_per_elem;
         elem_vert < (elem1 + 1) * nverts_per_elem;
         ++elem_vert) {
      auto vert = elems2verts[elem_vert];
      for (auto vert_elem = verts2elems.a2ab[vert];
           vert_elem < verts2elems.a2ab[vert + 1];
           ++vert_elem) {
        auto elem2 = verts2elems.ab2b[vert_elem];
        n = add_unique(stack, n, elem2, stack_max);
      }
    }
  }
  return n;
}

DEVICE Int find_distance_3_keys(
      Graph const& keys2buf_elems,
      LOs const& elems2verts,
      Int nverts_per_elem,
      Graph const& verts2elems,
      LOs const& elems2keys,
      Int nkeys_per_elem,
      Read<I8> const& indset,
      LO key,
      LO* stack,
      Int stack_max) {
  if (!indset[key]) return 0;
  Int n = 0;
  for (auto key_buf_elem = keys2buf_elems.a2ab[key];
       key_buf_elem < keys2buf_elems.a2ab[key + 1];
       ++key_buf_elem) {
    auto elem1 = keys2buf_elems.ab2b[key_buf_elem];
    for (auto elem_vert = elem1 * nverts_per_elem;
         elem_vert < (elem1 + 1) * nverts_per_elem;
         ++elem_vert) {
      auto vert = elems2verts[elem_vert];
      for (auto vert_elem = verts2elems.a2ab[vert];
           vert_elem < verts2elems.a2ab[vert + 1];
           ++vert_elem) {
        auto elem2 = verts2elems.ab2b[vert_elem];
        for (auto elem_key = elem2 * nkeys_per_elem;
             elem_key < (elem2 + 1) * nkeys_per_elem;
             ++elem_key) {
          auto key2 = elems2keys[elem_key];
          if (key2 != key && indset[key2]) {
            n = add_unique(stack, n, key2, stack_max);
          }
        }
      }
    }
  }
  return n;
}

} // end anonymous namespace

Graph get_buffered_elems(Mesh* mesh, Int key_dim,
    Read<I8> unbuffered_indset) {
  auto keys2elems = mesh->ask_up(key_dim, mesh->dim());
  auto elems2verts = mesh->ask_elem_verts();
  auto nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto nkeys = mesh->nents(key_dim);
  auto degrees_w = Write<LO>(nkeys);
  auto count = LAMBDA(LO key) {
    constexpr Int stack_max = 256;
    LO stack[stack_max];
    degrees_w[key] = find_distance_2_elems(keys2elems, elems2verts,
        nverts_per_elem, verts2elems,
        unbuffered_indset, key, stack, stack_max);
  };
  parallel_for(nkeys, count);
  auto degrees = LOs(degrees_w);
  auto offsets = offset_scan(degrees);
  auto edges = Write<LO>(offsets.last());
  auto fill = LAMBDA(LO key) {
    auto begin = offsets[key];
    auto end = offsets[key + 1];
    auto stack_max = end - begin;
    auto stack = edges.data() + begin;
    find_distance_2_elems(keys2elems, elems2verts,
        nverts_per_elem, verts2elems,
        unbuffered_indset, key, stack, stack_max);
  };
  parallel_for(nkeys, fill);
  return Graph(offsets, edges);
}

Graph get_buffered_conflicts(Mesh* mesh, Int key_dim,
    Read<I8> unbuffered_indset) {
  auto keys2buf_elems = get_buffered_elems(mesh, key_dim, unbuffered_indset);
  auto elems2verts = mesh->ask_elem_verts();
  auto nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto elems2keys = mesh->ask_down(mesh->dim(), key_dim).ab2b;
  auto nkeys_per_elem = simplex_degrees[mesh->dim()][key_dim];
  auto nkeys = mesh->nents(key_dim);
  auto degrees_w = Write<LO>(nkeys);
  auto count = LAMBDA(LO key) {
    constexpr Int stack_max = 32;
    LO stack[stack_max];
    degrees_w[key] = find_distance_3_keys(keys2buf_elems, elems2verts,
        nverts_per_elem, verts2elems, elems2keys, nkeys_per_elem,
        unbuffered_indset, key, stack, stack_max);
  };
  parallel_for(nkeys, count);
  auto degrees = LOs(degrees_w);
  auto offsets = offset_scan(degrees);
  auto edges = Write<LO>(offsets.last());
  auto fill = LAMBDA(LO key) {
    auto begin = offsets[key];
    auto end = offsets[key + 1];
    auto stack_max = end - begin;
    auto stack = edges.data() + begin;
    find_distance_3_keys(keys2buf_elems, elems2verts,
        nverts_per_elem, verts2elems, elems2keys, nkeys_per_elem,
        unbuffered_indset, key, stack, stack_max);
  };
  parallel_for(nkeys, fill);
  return Graph(offsets, edges);
}

Read<I8> find_buffered_indset(
    Mesh* mesh, Int key_dim,
    Reals qualities,
    Read<I8> unbuffered_indset) {
  auto conflicts = get_buffered_conflicts(mesh, key_dim, unbuffered_indset);
  return find_indset(mesh, key_dim, conflicts, qualities, unbuffered_indset);
}

} // end namespace Omega_h
