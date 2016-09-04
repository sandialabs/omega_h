#include "internal.hpp"
#include "loop.hpp"
#include "scan.hpp"

namespace Omega_h {

namespace {

DEVICE Int add_unique(LO* stack, Int n, LO v, Int stack_max) {
  for (Int i = 0; i < n; ++i)
    if (stack[i] == v) return n;
  CHECK(n < stack_max);
  stack[n++] = v;
  return n;
}

DEVICE Int distance_3_from_node(Graph const& g,
    Read<I8> const& indset, LO v, LO* stack, Int stack_max) {
  if (!indset[v]) return 0;
  Int n = 0;
  for (auto ab1 = g.a2ab[v]; ab1 < g.a2ab[v + 1]; ++ab1) {
    auto b1 = g.ab2b[ab1];
    for (auto ab2 = g.a2ab[b1]; ab2 < g.a2ab[b1 + 1]; ++ab2) {
      auto b2 = g.ab2b[ab2];
      for (auto ab3 = g.a2ab[b2]; ab3 < g.a2ab[b2 + 1]; ++ab3) {
        auto b3 = g.ab2b[ab3];
        if (v != b3 && indset[b3]) {
          n = add_unique(stack, n, b3, stack_max);
        }
      }
    }
  }
  return n;
}

LOs count(Graph g, Read<I8> indset) {
  auto nnodes = g.a2ab.size() - 1;
  auto degrees_w = Write<LO>(nnodes);
  auto f = LAMBDA(LO v) {
    constexpr Int max_adj = 32;
    LO stack[max_adj];
    degrees_w[v] = distance_3_from_node(g, indset, v, stack, max_adj);
  };
  parallel_for(nnodes, f);
  return degrees_w;
}

LOs fill(Graph g, Read<I8> indset, LOs offsets) {
  auto nnodes = g.a2ab.size() - 1;
  auto out = Write<LO>(offsets.last());
  auto f = LAMBDA(LO v) {
    auto begin = offsets[v];
    auto end = offsets[v + 1];
    auto stack_max = end - begin;
    distance_3_from_node(g, indset, v, out.data() + begin, stack_max);
  };
  parallel_for(nnodes, f);
  return out;
}

} // end anonymous namespace

Graph get_buffered_conflict_graph(Graph unbuffered_conflicts,
    Read<I8> unbuffered_indset) {
  auto g = unbuffered_conflicts;
  auto indset = unbuffered_indset;
  auto degrees = count(g, indset);
  auto a2ab = offset_scan(degrees);
  auto ab2b = fill(g, indset, a2ab);
  return Graph{a2ab, ab2b};
}

} // end namespace Omega_h
