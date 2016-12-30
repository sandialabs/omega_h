#ifndef OMEGA_H_MOTION_CHOICES_HPP
#define OMEGA_H_MOTION_CHOICES_HPP

#include "Omega_h.hpp"

namespace Omega_h {

struct MotionChoices {
  Reals qualities;
  Reals coordinates;
};

MotionChoices get_motion_choices(Mesh* mesh, AdaptOpts const& opts,
    LOs cands2verts);

}

#endif
