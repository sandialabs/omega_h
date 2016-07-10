#ifndef COLLAPSE_HPP
#define COLLAPSE_HPP

#include "internal.hpp"

namespace osh {

enum {
  DONT_COLLAPSE = 0x0,
  COLLAPSE_LEFT = 0x1,
  COLLAPSE_RIGHT = 0x2,
  COLLAPSE_BOTH = 0x3
};

INLINE I8 collapses(I8 code, Int col_v) { return (code & (1 << col_v)) != 0; }

INLINE I8 do_collapse(I8 code, Int col_v) {
  return static_cast<I8>(code | (1 << col_v));
}

INLINE I8 dont_collapse(I8 code, Int col_v) {
  return static_cast<I8>(code & ~(1 << col_v));
}

Read<I8> check_collapse_class(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);

Read<I8> check_collapse_exposure(Mesh* mesh, LOs cands2edges,
                                 Read<I8> cand_codes);

}  // end namespace osh

#endif
