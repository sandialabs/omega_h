#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "internal.hpp"

namespace osh {

namespace surf {

Reals get_side_normals(Mesh* mesh, LOs surf_side2side);
Reals get_hinge_angles(Mesh* mesh, Reals surf_side_normals,
                       LOs surf_hinge2hinge, LOs side2surf_side);
Reals get_vert_normals(Mesh* mesh, LOs surf_side2side, Reals surf_side_normals,
                       LOs surf_vert2vert);

}  // end namespace surf

}  // end namespace osh

#endif
