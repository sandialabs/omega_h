#ifndef OMEGA_H_SURFACE_HPP
#define OMEGA_H_SURFACE_HPP

#include "internal.hpp"

namespace Omega_h {

Reals get_side_normals(Mesh* mesh, LOs surf_side2side);
Reals get_curv_edge_tangents(Mesh* mesh, LOs curv_edge2edge);
Reals get_hinge_angles(Mesh* mesh, Reals surf_side_normals,
    LOs surf_hinge2hinge, LOs side2surf_side);
Reals get_side_vert_normals(Mesh* mesh, LOs surf_side2side,
    Reals surf_side_normals, LOs surf_vert2vert);
Reals get_curv_vert_tangents(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_tangents, LOs curv_vert2vert);
Reals get_surf_tri_IIs(Mesh* mesh, LOs surf_tris2tri, Reals surf_tri_normals,
    LOs surf_verts2vert, Reals surf_vert_normals);
Reals get_surf_vert_IIs(Mesh* mesh, LOs surf_tris2tri, Reals surf_tri_normals,
    Reals surf_tri_IIs, LOs surf_verts2vert, Reals surf_vert_normals);
Reals get_curv_edge_curvatures(Mesh* mesh, LOs curv_edges2edge,
    Reals curv_edge_tangents, LOs curv_verts2vert, Reals curv_vert_tangents);
Reals get_curv_vert_curvatures(Mesh* mesh, LOs curv_edges2edge,
    Reals curv_edge_curvatures, LOs curv_verts2vert);
Reals get_corner_vert_curvatures(Mesh* mesh, Write<Real> vert_curvatures_w);
Reals get_vert_curvatures(Mesh* mesh);

}  // end namespace Omega_h

#endif
