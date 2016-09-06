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
void transfer_conserve_r3d(Mesh* old_mesh, Mesh* new_mesh, Int key_dim,
    LOs keys2kds, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents);
void transfer_conserve_r3d_refine(Mesh* old_mesh, Mesh* new_mesh,
    LOs keys2edges, LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
    LOs same_ents2new_ents);
bool needs_buffer_layers(Mesh* mesh);
Graph get_buffered_conflict_graph(Mesh* mesh, Int key_dim,
    Read<I8> unbuffered_indset);
Graph get_buffered_elems(Mesh* mesh, Int key_dim,
    Read<I8> buffered_indset);
Read<I8> find_buffered_indset(
    Mesh* mesh, Int key_dim,
    Reals qualities,
    Read<I8> unbuffered_indset);
}

#endif
