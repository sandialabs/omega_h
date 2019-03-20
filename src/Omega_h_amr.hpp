#ifndef OMEGA_H_AMR_HPP
#define OMEGA_H_AMR_HPP

#include <Omega_h_adapt.hpp>
#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

namespace amr {

OMEGA_H_INLINE constexpr Int code_which_child(I8 code) { return (code >> 2); }

OMEGA_H_INLINE constexpr Int code_parent_dim(I8 code) { return code & 3; }

OMEGA_H_INLINE constexpr I8 make_code(Int which_child, Int parent_dim) {
  return static_cast<I8>((which_child << 2) | parent_dim);
}

void remove_non_leaf_uses(Mesh* mesh);
Bytes enforce_2to1_refine(Mesh* mesh, Int bridge_dim, Bytes elems_are_marked);
void refine(Mesh* mesh, Bytes elems_are_marked, TransferOpts xfer_opts);
void derefine(Mesh* mesh, Bytes elems_are_marked, TransferOpts xfer_opts);

}  // namespace amr

}  // namespace Omega_h

#endif
