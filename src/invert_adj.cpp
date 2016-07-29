#include "adjacency.hpp"

#include "algorithm.hpp"
#include "align.hpp"
#include "loop.hpp"
#include "map.hpp"

namespace osh {

static void order_by_globals(LOs l2lh, Write<LO> lh2h, Write<I8> codes,
                             Read<GO> hg) {
  LO nl = l2lh.size() - 1;
  auto f = LAMBDA(LO l) {
    LO begin = l2lh[l];
    LO end = l2lh[l + 1];
    for (auto i = begin; i < end; ++i) {
      auto k = i;
      auto i_key = hg[lh2h[i]];
      for (auto j = i + 1; j < end; ++j) {
        auto j_key = hg[lh2h[j]];
        if (j_key < i_key) {
          k = j;
          i_key = j_key;
        }
      }
      auto tmp_h = lh2h[i];
      auto tmp_code = codes[i];
      lh2h[i] = lh2h[k];
      codes[i] = codes[k];
      lh2h[k] = tmp_h;
      codes[k] = tmp_code;
    }
  };
  parallel_for(nl, f);
}

Adj invert_adj(Adj down, Int nlows_per_high, LO nlows, Read<GO> high_globals) {
  auto l2hl = invert_map_by_atomics(down.ab2b, nlows);
  auto l2lh = l2hl.a2ab;
  auto lh2hl = l2hl.ab2b;
  LO nlh = lh2hl.size();
  Read<I8> down_codes(down.codes);
  Write<LO> lh2h(nlh);
  Write<I8> codes(nlh);
  if (down_codes.exists()) {
    auto f = LAMBDA(LO lh) {
      LO hl = lh2hl[lh];
      LO h = hl / nlows_per_high;
      lh2h[lh] = h;
      Int which_down = hl % nlows_per_high;
      auto down_code = down_codes[hl];
      bool is_flipped = code_is_flipped(down_code);
      Int rotation = code_rotation(down_code);
      codes[lh] = make_code(is_flipped, rotation, which_down);
    };
    parallel_for(nlh, f);
  } else {
    auto f = LAMBDA(LO lh) {
      LO hl = lh2hl[lh];
      LO h = hl / nlows_per_high;
      lh2h[lh] = h;
      Int which_down = hl % nlows_per_high;
      codes[lh] = make_code(0, 0, which_down);
    };
    parallel_for(nlh, f);
  }
  order_by_globals(l2lh, lh2h, codes, high_globals);
  return Adj(l2lh, lh2h, codes);
}

}  // end namespace osh
