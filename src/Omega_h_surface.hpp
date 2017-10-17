#ifndef OMEGA_H_SURFACE_HPP
#define OMEGA_H_SURFACE_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

Reals get_side_vectors(Mesh* mesh, LOs surf_side2side);
Reals get_curv_edge_tangents(Mesh* mesh, LOs curv_edge2edge);
Reals get_hinge_angles(Mesh* mesh, Reals surf_side_normals,
    LOs surf_hinge2hinge, LOs side2surf_side);
Reals get_side_vert_normals(Mesh* mesh, LOs surf_side2side,
    Reals surf_side_normals, LOs surf_vert2vert);
Reals get_curv_vert_tangents(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_tangents, LOs curv_vert2vert);
Reals get_surf_tri_IIs(Mesh* mesh, LOs surf_tri2tri, Reals surf_tri_normals,
    LOs surf_vert2vert, Reals surf_vert_normals);
Reals get_surf_vert_IIs(Mesh* mesh, LOs surf_tri2tri, Reals surf_tri_normals,
    Reals surf_tri_IIs, LOs surf_vert2vert, Reals surf_vert_normals);
Reals get_curv_edge_curvatures(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_tangents, LOs curv_vert2vert, Reals curv_vert_tangents);
Reals get_curv_vert_curvatures(Mesh* mesh, LOs curv_edge2edge,
    Reals curv_edge_curvatures, LOs curv_vert2vert);

struct SurfaceInfo {
  LOs surf_vert2vert;
  Reals surf_vert_normals;
  Reals surf_vert_IIs;
  LOs curv_vert2vert;
  Reals curv_vert_tangents;
  Reals curv_vert_curvatures;
};

SurfaceInfo get_surface_info(Mesh* mesh);

Reals get_vert_curvatures(Mesh* mesh, SurfaceInfo surface_info);

}  // end namespace Omega_h

#endif
