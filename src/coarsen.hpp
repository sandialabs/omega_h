#ifndef COARSEN_HPP
#define COARSEN_HPP

#include "host_few.hpp"

namespace Omega_h {

/* This is (1/sqrt(2)), but sqrt() is not constexpr */
constexpr Real min_length_desired = 0.7071067811865475;

HostFew<Read<I8>, 4> mark_dead_ents(
    Mesh* mesh, LOs rails2edges, Read<I8> rail_col_dirs);

Adj find_coarsen_domains(
    Mesh* mesh, LOs keys2verts, Int ent_dim, Read<I8> ents_are_dead);

Reals coarsen_qualities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);

Read<I8> filter_coarsen_min_qual(
    Read<I8> cand_codes, Reals cand_quals, Real min_qual);

Read<I8> filter_coarsen_improve(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes, Reals cand_quals);

Read<I8> prevent_overshoot(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);

void choose_rails(Mesh* mesh, LOs cands2edges,
    Read<I8> cand_edge_codes, Reals cand_edge_quals, Read<I8>* verts_are_cands,
    Reals* vert_quals, Read<GO>* vert_rails);
void find_rails(Mesh* mesh, LOs keys2verts, Read<GO> verts2rail,
    LOs* rails2edges, Read<I8>* rail_col_dirs);

LOs get_verts_onto(Mesh* mesh, LOs rails2edges, Read<I8> rail_col_dirs);

LOs coarsen_topology(Mesh* mesh, LOs keys2verts_onto, Int dom_dim,
    Adj keys2doms, LOs old_verts2new_verts);

bool coarsen_by_size(Mesh* mesh, AdaptOpts const& opts);

bool coarsen_slivers(Mesh* mesh, AdaptOpts const& opts);

}  // end namespace Omega_h

#endif
