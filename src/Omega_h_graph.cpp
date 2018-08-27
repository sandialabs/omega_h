#include "Omega_h_graph.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_map.hpp"

namespace Omega_h {

LO Graph::nnodes() const { return a2ab.size() - 1; }

LO Graph::nedges() const { return ab2b.size(); }

Graph add_edges(Graph g1, Graph g2) {
  OMEGA_H_TIME_FUNCTION;
  auto v2e1 = g1.a2ab;
  auto e2v1 = g1.ab2b;
  auto v2e2 = g2.a2ab;
  auto e2v2 = g2.ab2b;
  auto nv = v2e1.size() - 1;
  auto deg1 = get_degrees(v2e1, "deg1");
  auto deg2 = get_degrees(v2e2, "deg2");
  auto deg = add_each(deg1, deg2, "deg1+2");
  auto v2e = offset_scan(deg, v2e1.name());
  Write<LO> e2v(v2e.last(), e2v1.name());
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto begin1 = v2e1[v];
    auto end1 = v2e1[v + 1];
    auto begin2 = v2e2[v];
    auto end2 = v2e2[v + 1];
    auto begin = v2e[v];
    auto end = v2e[v + 1];
    auto k = begin;
    for (auto j = begin1; j < end1; ++j) e2v[k++] = e2v1[j];
    for (auto j = begin2; j < end2; ++j) e2v[k++] = e2v2[j];
    OMEGA_H_CHECK(k == end);
  };
  parallel_for(nv, f, "add_edges");
  return Graph(v2e, e2v);
}

Graph unmap_graph(LOs a2b, Graph b2c) {
  auto b2bc = b2c.a2ab;
  auto bc2c = b2c.ab2b;
  auto b_degrees = get_degrees(b2bc);
  auto a_degrees = read(unmap(a2b, b_degrees, 1));
  auto a2ac = offset_scan(a_degrees);
  auto na = a2b.size();
  Write<LO> ac2c(a2ac.last());
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto b = a2b[a];
    auto bc = b2bc[b];
    for (auto ac = a2ac[a]; ac < a2ac[a + 1]; ++ac) {
      ac2c[ac] = bc2c[bc++];
    }
  };
  parallel_for(na, f, "unmap_graph");
  return Graph(a2ac, ac2c);
}

template <typename T>
Read<T> graph_reduce(Graph a2b, Read<T> b_data, Int width, Omega_h_Op op) {
  auto a2ab = a2b.a2ab;
  auto ab2b = a2b.ab2b;
  auto ab_data = read(unmap(ab2b, b_data, width));
  return fan_reduce(a2ab, ab_data, width, op);
}

Reals graph_weighted_average_arc_data(
    Graph a2b, Reals ab_weights, Reals ab_data, Int width) {
  auto a2ab = a2b.a2ab;
  auto ab2b = a2b.ab2b;
  auto nab = a2ab.last();
  OMEGA_H_CHECK(ab_weights.size() == nab);
  OMEGA_H_CHECK(ab_data.size() % width == 0);
  auto total_weights = fan_reduce(a2ab, ab_weights, 1, OMEGA_H_SUM);
  Reals weighted_ab_data = multiply_each(ab_data, ab_weights);
  auto weighted_sums = fan_reduce(a2ab, weighted_ab_data, width, OMEGA_H_SUM);
  return divide_each(weighted_sums, total_weights);
}

Reals graph_weighted_average(
    Graph a2b, Reals ab_weights, Reals b_data, Int width) {
  auto ab2b = a2b.ab2b;
  auto ab_data = unmap(ab2b, b_data, width);
  return graph_weighted_average_arc_data(a2b, ab_weights, ab_data, width);
}

struct FilteredGraph {
  Graph g;
  LOs kept2old;
};

Graph filter_graph(Graph g, Read<I8> keep_edge) {
  auto degrees = fan_reduce(g.a2ab, keep_edge, 1, OMEGA_H_SUM);
  auto offsets = offset_scan(degrees);
  auto kept2old = collect_marked(keep_edge);
  auto edges = unmap(kept2old, g.ab2b, 1);
  return Graph(offsets, edges);
}

bool operator==(Graph a, Graph b) {
  return a.a2ab == b.a2ab && a.ab2b == b.ab2b;
}

Graph identity_graph(LO nnodes) {
  return {LOs(nnodes + 1, 0, 1), LOs(nnodes, 0, 1)};
}

Graph add_self_edges(Graph g) {
  return add_edges(identity_graph(g.nnodes()), g);
}

template <typename T>
void map_into(Read<T> a_data, Graph a2b, Write<T> b_data, Int width) {
  auto ab_data = expand(a_data, a2b.a2ab, width);
  map_into(ab_data, a2b.ab2b, b_data, width);
}

template <typename T>
Read<T> map_onto(Read<T> a_data, Graph a2b, LO nb, T init_val, Int width) {
  auto out = Write<T>(nb * width, init_val);
  map_into(a_data, a2b, out, width);
  return out;
}

#define INST(T)                                                                \
  template Read<T> graph_reduce(Graph, Read<T>, Int, Omega_h_Op);              \
  template void map_into(                                                      \
      Read<T> a_data, Graph a2b, Write<T> b_data, Int width);                  \
  template Read<T> map_onto(Read<T> a_data, Graph a2b, LO nb, T, Int width);
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
