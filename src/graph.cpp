#include "graph.hpp"

#include "array.hpp"
#include "loop.hpp"
#include "map.hpp"
#include "scan.hpp"

namespace osh {

Graph add_edges(Graph g1, Graph g2) {
  auto v2e1 = g1.a2ab;
  auto e2v1 = g1.ab2b;
  auto v2e2 = g2.a2ab;
  auto e2v2 = g2.ab2b;
  auto nv = v2e1.size() - 1;
  auto deg1 = get_degrees(v2e1);
  auto deg2 = get_degrees(v2e2);
  auto deg = add_each(deg1, deg2);
  auto v2e = offset_scan(deg);
  Write<LO> e2v(v2e.last());
  auto f = LAMBDA(LO v) {
    auto begin1 = v2e1[v];
    auto end1 = v2e1[v + 1];
    auto begin2 = v2e2[v];
    auto end2 = v2e2[v + 1];
    auto begin = v2e[v];
    auto end = v2e[v + 1];
    auto k = begin;
    for (auto j = begin1; j < end1; ++j) e2v[k++] = e2v1[j];
    for (auto j = begin2; j < end2; ++j) e2v[k++] = e2v2[j];
    CHECK(k == end);
  };
  parallel_for(nv, f);
  return Graph(v2e, e2v);
}

Graph unmap_graph(LOs a2b, Graph b2c) {
  auto b2bc = b2c.a2ab;
  auto bc2c = b2c.ab2b;
  auto b_degrees = get_degrees(b2bc);
  auto a_degrees = unmap(a2b, b_degrees, 1);
  auto a2ac = offset_scan(a_degrees);
  auto na = a2b.size();
  Write<LO> ac2c(a2ac.last());
  auto f = LAMBDA(LO a) {
    auto b = a2b[a];
    auto bc = b2bc[b];
    for (auto ac = a2ac[a]; ac < a2ac[a + 1]; ++ac) {
      ac2c[ac] = bc2c[bc++];
    }
  };
  parallel_for(na, f);
  return Graph(a2ac, ac2c);
}

Adj unmap_adjacency(LOs a2b, Adj b2c) {
  auto b2bc = b2c.a2ab;
  auto bc2c = b2c.ab2b;
  auto bc_codes = b2c.codes;
  auto b_degrees = get_degrees(b2bc);
  auto a_degrees = unmap(a2b, b_degrees, 1);
  auto a2ac = offset_scan(a_degrees);
  auto na = a2b.size();
  auto nac = a2ac.last();
  Write<LO> ac2c(nac);
  auto ac_codes = Write<I8>(nac);
  auto f = LAMBDA(LO a) {
    auto b = a2b[a];
    auto bc = b2bc[b];
    for (auto ac = a2ac[a]; ac < a2ac[a + 1]; ++ac) {
      ac2c[ac] = bc2c[bc];
      ac_codes[ac] = bc_codes[bc];
      ++bc;
    }
  };
  parallel_for(na, f);
  return Adj(a2ac, ac2c, ac_codes);
}

template <typename T>
Read<T> graph_reduce(Graph a2b, Read<T> b_data, Int width, osh_op op) {
  auto a2ab = a2b.a2ab;
  auto ab2b = a2b.ab2b;
  auto ab_data = unmap(ab2b, b_data, width);
  return fan_reduce(a2ab, ab_data, width, op);
}

Reals graph_weighted_average_arc_data(
    Graph a2b, Reals ab_weights, Reals ab_data, Int width) {
  auto a2ab = a2b.a2ab;
  auto ab2b = a2b.ab2b;
  auto nab = a2ab.last();
  CHECK(ab_weights.size() == nab);
  CHECK(ab_data.size() % width == 0);
  auto total_weights = fan_reduce(a2ab, ab_weights, 1, OSH_SUM);
  auto weighted_ab_data = multiply_each(ab_data, ab_weights);
  auto weighted_sums = fan_reduce(a2ab, weighted_ab_data, width, OSH_SUM);
  return divide_each(weighted_sums, total_weights);
}

Reals graph_weighted_average(
    Graph a2b, Reals ab_weights, Reals b_data, Int width) {
  auto ab2b = a2b.ab2b;
  auto ab_data = unmap(ab2b, b_data, width);
  return graph_weighted_average_arc_data(a2b, ab_weights, ab_data, width);
}

#define INST(T) template Read<T> graph_reduce(Graph, Read<T>, Int, osh_op);
INST(I8)
INST(I32)
INST(Real)
#undef INST

}  // end namespace osh
