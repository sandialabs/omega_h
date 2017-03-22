#ifndef BOX_HPP
#define BOX_HPP

#include "internal.hpp"

namespace Omega_h {

void make_1d_box(
    Real x, LO nx, LOs* ev2v_out, Reals* coords_out);
void make_2d_box(
    Real x, Real y, LO nx, LO ny, LOs* qv2v_out, Reals* coords_out);
void make_3d_box(Real x, Real y, Real z, LO nx, LO ny, LO nz, LOs* hv2v_out,
    Reals* coords_out);
void set_box_class_ids(
    Mesh* mesh, Real x, Real y, Real z, Int nx, Int ny, Int nz);

}  // end namespace Omega_h

#endif
