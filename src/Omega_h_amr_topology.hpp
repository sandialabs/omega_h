#ifndef OMEGA_H_AMR_TOPOLOGY_HPP
#define OMEGA_H_AMR_TOPOLOGY_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_few.hpp>
#include <Omega_h_template_up.hpp>

namespace Omega_h {

class Mesh;

void mark_amr(Mesh* mesh, Bytes elems_are_marked);

Few<LO, 4> count_amr(Mesh* mesh);

LOs get_amr_topology(Mesh* mesh, Int child_dim, Int num_children,
    Few<LOs, 4> mods2mds, Few<LOs, 4> mds2mods, Few<LOs, 4> mods2midverts,
    LOs old_verts2new_verts);

}  // namespace Omega_h

#endif
