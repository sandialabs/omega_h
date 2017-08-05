#ifndef OMEGA_H_INDSET_HPP
#define OMEGA_H_INDSET_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_graph.hpp>
#include <Omega_h_dist.hpp>

namespace Omega_h {

class Mesh;

GOs find_indset(
    Graph graph,
    Int distance,
    Bytes candidates,
    Reals qualities,
    GOs globals, 
    Dist owners2copies);

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Graph graph, Reals qualities, Bytes candidates);
Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Reals qualities, Bytes candidates);

}  // end namespace Omega_h

#endif
