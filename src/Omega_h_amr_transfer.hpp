#ifndef OMEGA_H_AMR_TRANSFER_HPP
#define OMEGA_H_AMR_TRANSFER_HPP

#include <Omega_h_adapt.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_few.hpp>

namespace Omega_h {

void amr_transfer_linear_interp(Mesh* old_mesh, Mesh* new_mesh,
    Few<LOs, 4> mods2mds, Few<LOs, 4> mods2midverts, LOs same_ents2old_ents,
    LOs same_ents2new_ents, TransferOpts opts);

void amr_transfer_levels(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
    LOs same_ents2old_ents, LOs same_ents2new_ents, Few<LOs, 4> mods2mds,
    LOs prods2new_ents);

void amr_transfer_leaves(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
    LOs same_ents2old_ents, LOs same_ents2new_ents, Few<LOs, 4> mods2mds,
    LOs prods2new_ents, LOs old_ents2new_ents);

void amr_transfer_parents(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim,
    LOs same_ents2old_ents, LOs same_ents2new_ents, Few<LOs, 4> mods2mds,
    LOs prods2new_ents, LOs old_ents2new_ents);

}  // namespace Omega_h

#endif
