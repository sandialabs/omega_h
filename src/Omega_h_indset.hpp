#ifndef OMEGA_H_INDSET_HPP
#define OMEGA_H_INDSET_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_graph.hpp>

namespace Omega_h {

class Mesh;

Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Graph graph, Reals quality, Read<I8> candidates);
Read<I8> find_indset(
    Mesh* mesh, Int ent_dim, Reals quality, Read<I8> candidates);

}  // end namespace Omega_h

#endif
