#ifndef OMEGA_H_TRANSFER_CONSERVE_HPP
#define OMEGA_H_TRANSFER_CONSERVE_HPP

#include <Omega_h_adapt.hpp>
#include <Omega_h_adj.hpp>

namespace Omega_h {

class Mesh;

void transfer_conserve_refine(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
    LOs same_ents2old_ents, LOs same_ents2new_ents);
void transfer_densities_and_conserve_swap(Mesh* old_mesh,
    TransferOpts const& opts, Mesh* new_mesh, LOs keys2edges, LOs keys2prods,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents);
void transfer_densities_and_conserve_coarsen(Mesh* old_mesh,
    TransferOpts const& opts, Mesh* new_mesh, LOs keys2verts, Adj keys2doms,
    LOs prods2new_ents, LOs same_ents2old_ents, LOs same_ents2new_ents);
void transfer_conserve_motion(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, LOs keys2verts, Graph keys2elems, LOs same_ents2old_ents,
    LOs same_ents2new_ents);
void setup_conservation_tags(Mesh* mesh, AdaptOpts const& opts);
void correct_integral_errors(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
