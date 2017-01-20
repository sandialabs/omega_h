#ifndef OMEGA_H_MOTION_CHOICES_HPP
#define OMEGA_H_MOTION_CHOICES_HPP

#include "Omega_h.hpp"

namespace Omega_h {

struct MotionChoices {
  Read<I8> cands_did_move;
  Reals quals;
  Reals new_sol;
};

void unpack_linearized_fields(Mesh* old_mesh, Mesh* new_mesh, Reals data);

MotionChoices get_motion_choices(Mesh* mesh, AdaptOpts const& opts,
    LOs cands2verts);

}

#endif
