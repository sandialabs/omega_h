#ifndef OMEGA_H_MOTION_HPP
#define OMEGA_H_MOTION_HPP

#include "Omega_h.hpp"

namespace Omega_h {

struct LinearPack {
  Reals data;
  Int ncomps;
  Int metric_offset;
  Int coords_offset;
};

LinearPack pack_linearized_fields(Mesh* mesh, TransferOpts const& opts);
void unpack_linearized_fields(Mesh* old_mesh, TransferOpts const& opts,
    Mesh* new_mesh, Reals data, Read<I8> verts_are_keys);

struct MotionChoices {
  Read<I8> did_move;
  Reals new_quals;
  Reals new_coords;
};

MotionChoices get_motion_choices(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2verts);

bool move_verts_for_quality(Mesh* mesh, AdaptOpts const& opts);

}  // namespace Omega_h

#endif
