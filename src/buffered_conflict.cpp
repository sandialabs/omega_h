#include "internal.hpp"
#include "loop.hpp"
#include "scan.hpp"
#include "simplices.hpp"

namespace Omega_h {

namespace {

DEVICE Int add_unique(LO* stack, Int n, LO e, Int stack_max) {
  for (Int i = 0; i < n; ++i)
    if (stack[i] == e) return n;
  CHECK(n < stack_max);
  stack[n++] = e;
  return n;
}

struct KeysDistance3 {
  enum { max_adj = 32 };
  Graph keys2elems;
  LOs elems2verts;
  Int nverts_per_elem;
  Graph verts2elems;
  LOs elems2keys;
  Int nkeys_per_elem;
  KeysDistance3(Mesh* mesh, Int key_dim) {
    keys2elems = mesh->ask_up(key_dim, mesh->dim());
    elems2verts = mesh->ask_elem_verts();
    nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
    verts2elems = mesh->ask_up(VERT, mesh->dim());
    elems2keys = mesh->ask_down(mesh->dim(), key_dim).ab2b;
    nkeys_per_elem = simplex_degrees[mesh->dim()][key_dim];
  }
  DEVICE Int find(
      Read<I8> const& indset,
      LO key,
      LO* stack,
      Int stack_max) const {
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
          for (auto elem_key = elem2 * nkeys_per_elem;
               elem_key < (elem2 + 1) * nkeys_per_elem;
               ++elem_key) {
            auto key2 = elems2keys[elem_key];
            if (indset[key2]) {
              n = add_unique(stack, n, key2, stack_max);
            }
          }
        }
      }
    }
    return n;
  }
};

struct ElemsDistance2 {
  enum { max_adj = 256 };
  Graph keys2elems;
  LOs elems2verts;
  Int nverts_per_elem;
  Graph verts2elems;
  ElemsDistance2(Mesh* mesh, Int key_dim) {
    keys2elems = mesh->ask_up(key_dim, mesh->dim());
    elems2verts = mesh->ask_elem_verts();
    nverts_per_elem = simplex_degrees[mesh->dim()][VERT];
    verts2elems = mesh->ask_up(VERT, mesh->dim());
  }
  DEVICE Int find(
      Read<I8> const& indset,
      LO key,
      LO* stack,
      Int stack_max) const {
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
};

template <typename Op>
LOs count(Mesh* mesh, Int key_dim, Read<I8> indset) {
  auto nkeys = mesh->nents(key_dim);
  auto degrees_w = Write<LO>(nkeys);
  Op op(mesh, key_dim);
  auto f = LAMBDA(LO key) {
    LO stack[Op::max_adj];
    degrees_w[key] = op.find(indset, key, stack, Op::max_adj);
  };
  parallel_for(nkeys, f);
  return degrees_w;
}

template <typename Op>
LOs fill(Mesh* mesh, Int key_dim, Read<I8> indset, LOs offsets) {
  auto nkeys = mesh->nents(key_dim);
  auto out = Write<LO>(offsets.last());
  Op op(mesh, key_dim);
  auto f = LAMBDA(LO key) {
    auto begin = offsets[key];
    auto end = offsets[key + 1];
    auto stack_max = end - begin;
    auto stack = out.data() + begin;
    op.find(indset, key, stack, stack_max);
  };
  parallel_for(nkeys, f);
  return out;
}

template <typename Op>
Graph get_graph(Mesh* mesh, Int key_dim, Read<I8> indset) {
  auto degrees = count<KeysDistance3>(mesh, key_dim, indset);
  auto a2ab = offset_scan(degrees);
  auto ab2b = fill<KeysDistance3>(mesh, key_dim, indset, a2ab);
  return Graph{a2ab, ab2b};
}

} // end anonymous namespace

Graph get_buffered_conflict_graph(Mesh* mesh, Int key_dim,
    Read<I8> unbuffered_indset) {
  return get_graph<KeysDistance3>(mesh, key_dim, unbuffered_indset);
}

Graph get_buffered_elems(Mesh* mesh, Int key_dim,
    Read<I8> buffered_indset) {
  return get_graph<ElemsDistance2>(mesh, key_dim, buffered_indset);
}

} // end namespace Omega_h
