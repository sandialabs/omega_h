#ifndef TRANSFER_CONSERVE_HPP
#define TRANSFER_CONSERVE_HPP

#include "Omega_h.hpp"

namespace Omega_h {

void transfer_conserve_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
    LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents);
void transfer_conserve(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents);
void transfer_momentum_velocity(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_elems, LOs same_verts2old_verts,
    LOs same_verts2new_verts);

bool needs_buffer_layers(Mesh* mesh);

Graph get_buffered_elems(Mesh* mesh, Int key_dim, Read<I8> buffered_indset);
Graph get_buffered_conflicts(
    Mesh* mesh, Int key_dim, Graph kds2buf_elems, Read<I8> unbuffered_indset);

Read<I8> find_buffered_indset(
    Mesh* mesh, Int key_dim, Reals qualities, Read<I8> unbuffered_indset);

Graph get_closure_verts(Mesh* mesh, Graph keys2elems);

bool has_fixed_momentum_velocity(Mesh* mesh);
Read<I8> filter_coarsen_momentum_velocity(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);
Read<I8> filter_swap_momentum_velocity(Mesh* mesh, LOs cands2edges);

void do_momentum_velocity_elem_target(Mesh* donor_mesh, Mesh* target_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_elems);
void do_momentum_velocity_ghosted_target(Mesh* mesh);
}

#endif
