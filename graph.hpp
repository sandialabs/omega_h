struct Graph {
  Graph() {}
  Graph(LOs ab2b_):ab2b(ab2b_) {}
  Graph(LOs a2ab_, LOs ab2b_):a2ab(a2ab_),ab2b(ab2b_) {}
  LOs a2ab;
  LOs ab2b;
};

Graph add_edges(Graph g1, Graph g2);
Graph unmap_graph(LOs a2b, Graph b2c);
