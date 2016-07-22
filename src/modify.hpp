#ifndef MODIFY_HPP
#define MODIFY_HPP

#include "internal.hpp"

namespace osh {

LOs get_edge2rep_order(Mesh* mesh, Read<I8> edges_are_keys);

void modify_ents(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim, Int key_dim,
                 LOs keys2kds, LOs keys2prods, LOs prod_verts2verts,
                 LOs old_lows2new_lows, LOs* p_prods2new_ents,
                 LOs* p_same_ents2old_ents, LOs* p_same_ents2new_ents,
                 LOs* p_old_ents2new_ents);

void set_owners_by_indset(Mesh* mesh, Int key_dim, LOs keys2kds);

}  // end namespace osh

#endif
