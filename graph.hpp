Graph add_edges(Graph g1, Graph g2);
Graph unmap_graph(LOs a2b, Graph b2c);
Adj unmap_adjacency(LOs a2b, Adj b2c);
template <typename T>
Read<T> graph_reduce(Graph a2b, Read<T> b_data, Int width, osh_op op);
Reals graph_weighted_average(Graph a2b, Reals ab_weights, Reals b_data,
                             Int width);
