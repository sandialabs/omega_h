#ifndef OMEGA_H_MODIFY_HPP
#define OMEGA_H_MODIFY_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_few.hpp>
#include <Omega_h_graph.hpp>

namespace Omega_h {

class Mesh;

LOs get_rep2md_order_adapt(
    Mesh* mesh, Int key_dim, Int rep_dim, Bytes kds_are_keys);

Few<LOs, 4> get_rep2md_order(Mesh* mesh, Int rep_dim, Few<LOs, 4> mods2mds,
    Few<LOs, 4> mods2nprods, Few<bool, 4> mods_have_prods);

void modify_ents_adapt(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prod_verts2verts, LOs old_lows2new_lows,
    LOs* p_prods2new_ents, LOs* p_same_ents2old_ents, LOs* p_same_ents2new_ents,
    LOs* p_old_ents2new_ents);

void modify_ents(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    Few<LOs, 4> mods2mds, Few<Bytes, 4> mds_are_mods, Few<LOs, 4> mods2prods,
    LOs prod_verts2verts, LOs old_lows2new_lows, bool keep_mods,
    bool mods_can_be_shared, LOs* p_prods2new_ents, LOs* p_same_ents2old_ents,
    LOs* p_same_ents2new_ents, LOs* p_old_ents2new_ents);

void set_owners_by_indset(
    Mesh* mesh, Int key_dim, LOs keys2kds, Graph kds2elems);

}  // end namespace Omega_h

#endif
