#ifndef OMEGA_H_BUILD_HPP
#define OMEGA_H_BUILD_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>

namespace Omega_h {

class Mesh;
struct Remotes;

void build_verts_from_globals(Mesh* mesh, GOs vert_globals);
// builds only non-vertex entities, expects build_verts_from_globals to be
// called first
void build_ents_from_elems2verts(
    Mesh* mesh, LOs ev2v, GOs vert_globals, GOs elem_globals = GOs());
void build_from_elems2verts(Mesh* mesh, CommPtr comm, Omega_h_Family, Int edim,
    LOs ev2v, Read<GO> vert_globals);
void build_from_elems2verts(
    Mesh* mesh, Omega_h_Family family, Int edim, LOs ev2v, LO nverts);
void build_from_elems_and_coords(
    Mesh* mesh, Omega_h_Family family, Int edim, LOs ev2v, Reals coords);
Mesh build_box(CommPtr comm, Omega_h_Family family, Real x, Real y, Real z,
    LO nx, LO ny, LO nz, bool symmetric = false);
void build_box_internal(Mesh* mesh, Omega_h_Family family, Real x, Real y,
    Real z, LO nx, LO ny, LO nz, bool symmetric = false);

void add_ents2verts(
    Mesh* mesh, Int edim, LOs ev2v, GOs vert_globals, GOs elem_globals = GOs());
void resolve_derived_copies(CommPtr comm, Read<GO> verts2globs, Int deg,
    LOs* p_ent_verts2verts, Remotes* p_ents2owners);

void suggest_slices(
    GO total, I32 comm_size, I32 comm_rank, GO* p_begin, GO* p_end);
void assemble_slices(CommPtr comm, Omega_h_Family family, Int dim,
    GO global_nelems, GO elem_offset, GOs conn_in, GO global_nverts,
    GO vert_offset, Reals vert_coords, Dist* p_slice_elems2elems, LOs* conn_out,
    Dist* p_slice_verts2verts);

}  // end namespace Omega_h

#endif
