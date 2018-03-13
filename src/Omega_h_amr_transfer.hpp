#ifndef OMEGA_H_AMR_TRANSFER_HPP
#define OMEGA_H_AMR_TRANSFER_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_few.hpp>
#include <Omega_h_adapt.hpp>

namespace Omega_h {

void amr_transfer_linear_interp(Mesh* old_mesh, Mesh* new_mesh,
    Few<LOs, 4> mods2mds, Few<LOs, 4> mods2midverts,
    LOs same_ents2old_ents, LOs same_ents2new_ents,
    TransferOpts opts);

}

#endif
