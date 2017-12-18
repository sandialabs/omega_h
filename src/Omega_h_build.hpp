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
void build_from_elems2verts(
    Mesh* mesh, CommPtr comm, Int edim, LOs ev2v, Read<GO> vert_globals);
void build_from_elems2verts(Mesh* mesh, Int edim, LOs ev2v, LO nverts);
void build_from_elems_and_coords(Mesh* mesh, Int edim, LOs ev2v, Reals coords);
Mesh build_box(CommPtr comm, Real x, Real y, Real z, LO nx, LO ny, LO nz);
void build_box_internal(
    Mesh* mesh, Real x, Real y, Real z, LO nx, LO ny, LO nz);

void add_ents2verts(
    Mesh* mesh, Int edim, LOs ev2v, GOs vert_globals, GOs elem_globals = GOs());
void resolve_derived_copies(CommPtr comm, Read<GO> verts2globs, Int deg,
    LOs* p_ent_verts2verts, Remotes* p_ents2owners);

}  // end namespace Omega_h

#endif
