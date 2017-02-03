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

Graph get_closure_verts(Mesh* mesh, Graph keys2elems);

bool has_fixed_momentum_velocity(Mesh* mesh);
Read<I8> filter_coarsen_momentum_velocity(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);
Read<I8> filter_swap_momentum_velocity(Mesh* mesh, LOs cands2edges);

void do_momentum_velocity_part1(Mesh* donor_mesh, Mesh* target_mesh,
    Int key_dim, LOs keys2kds, LOs keys2prods, LOs prods2new_elems);
void do_momentum_velocity_part2(Mesh* mesh);

}  // end namespace Omega_h

#endif
