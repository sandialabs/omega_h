#ifndef TRANSFER_CONSERVE_HPP
#define TRANSFER_CONSERVE_HPP

#include "Omega_h.hpp"

namespace Omega_h {

void transfer_conserve_refine(Mesh* old_mesh, XferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents);
void transfer_conserve_swap(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents);
void transfer_conserve_coarsen(Mesh* old_mesh, XferOpts const& opts, Mesh* new_mesh,
    LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents);
void setup_conservation_tags(Mesh* mesh, AdaptOpts const& opts);
void correct_integral_errors(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
