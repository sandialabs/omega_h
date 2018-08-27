#include "Omega_h_parser_graph.hpp"

#include <iostream>

#include "Omega_h_std_vector.hpp"

namespace Omega_h {

ParserGraph make_graph_with_nnodes(int nnodes) {
  return ParserGraph(std::size_t(nnodes));
}

int get_nnodes(ParserGraph const& g) { return size(g); }

void add_edge(ParserGraph& g, int i, int j) { at(g, i).push_back(j); }

NodeEdges const& get_edges(ParserGraph const& g, int i) { return at(g, i); }

NodeEdges& get_edges(ParserGraph& g, int i) { return at(g, i); }

ParserGraph make_transpose(ParserGraph const& g) {
  auto nnodes = get_nnodes(g);
  auto transpose = make_graph_with_nnodes(nnodes);
  for (int i = 0; i < nnodes; ++i) {
    for (auto j : get_edges(g, i)) {
      add_edge(transpose, j, i);
    }
  }
  return transpose;
}

int at(ParserGraph const& g, int i, int j) { return at(at(g, i), j); }

std::ostream& operator<<(std::ostream& os, ParserGraph const& g) {
  for (int i = 0; i < get_nnodes(g); ++i) {
    os << i << ":";
    for (auto j : get_edges(g, i)) os << " " << j;
    os << '\n';
  }
  return os;
}

}  // namespace Omega_h
