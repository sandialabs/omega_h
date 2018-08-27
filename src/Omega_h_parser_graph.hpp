#ifndef OMEGA_H_GRAPH_HPP
#define OMEGA_H_GRAPH_HPP

#include <iosfwd>
#include <vector>

namespace Omega_h {

using NodeEdges = std::vector<int>;
using ParserGraph = std::vector<NodeEdges>;

ParserGraph make_graph_with_nnodes(int nnodes);
int get_nnodes(ParserGraph const& g);
void add_edge(ParserGraph& g, int i, int j);
NodeEdges const& get_edges(ParserGraph const& g, int i);
NodeEdges& get_edges(ParserGraph& g, int i);
ParserGraph make_transpose(ParserGraph const& g);
int at(ParserGraph const& g, int i, int j);
std::ostream& operator<<(std::ostream& os, ParserGraph const& g);

}  // namespace Omega_h

#endif
