#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "internal.hpp"

namespace osh {

Graph add_edges(Graph g1, Graph g2);
Graph unmap_graph(LOs a2b, Graph b2c);
Adj unmap_adjacency(LOs a2b, Adj b2c);
template <typename T>
Read<T> graph_reduce(Graph a2b, Read<T> b_data, Int width, osh_op op);
Reals graph_weighted_average_arc_data(Graph a2b, Reals ab_weights, Reals ab_data,
                             Int width);
Reals graph_weighted_average(Graph a2b, Reals ab_weights, Reals b_data,
                             Int width);

#define INST_DECL(T) \
  extern template Read<T> graph_reduce(Graph, Read<T>, Int, osh_op);
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(Real)
#undef INST_DECL

}  // end namespace osh

#endif
