#ifndef OMEGA_H_MODIFY_HPP
#define OMEGA_H_MODIFY_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_graph.hpp>

namespace Omega_h {

class Mesh;

LOs get_edge2rep_order(Mesh* mesh, Read<I8> edges_are_keys);

void modify_ents(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prod_verts2verts, LOs old_lows2new_lows,
    LOs* p_prods2new_ents, LOs* p_same_ents2old_ents, LOs* p_same_ents2new_ents,
    LOs* p_old_ents2new_ents);

void set_owners_by_indset(
    Mesh* mesh, Int key_dim, LOs keys2kds, Graph kds2elems);

}  // end namespace Omega_h

#endif
