#include "adjacency.hpp"

#include "algorithm.hpp"
#include "align.hpp"
#include "local_sort.hpp"
#include "loop.hpp"
#include "map.hpp"

#define MAX_UPWARD 100

namespace osh {

struct SortableAdj {
  struct value_type {
    LO high;
    I8 code;
  };
  typedef GO key_type;
  SortableAdj(Write<LO> const& lh2h, Write<I8> const& codes,
      Read<GO> const& hg, LO begin):hg_(hg) {
    lh2h_ptr_ = lh2h.data() + begin;
    codes_ptr_ = codes.data() + begin;
  }
  Read<GO> const& hg_;
  LO* lh2h_ptr_;
  I8* codes_ptr_;
  key_type key(Int i) const {
    return hg_[lh2h_ptr_[i]];
  }
  value_type value(Int i) const {
    return value_type{ lh2h_ptr_[i], codes_ptr_[i] };
  }
  void set(Int i, value_type const& x) const {
    lh2h_ptr_[i] = x.high;
    codes_ptr_[i] = x.code;
  }
};

static void order_by_globals(LOs l2lh, Write<LO> lh2h, Write<I8> codes,
                             Read<GO> hg) {
  LO nl = l2lh.size() - 1;
  auto f = LAMBDA(LO l) {
    LO begin = l2lh[l];
    LO end = l2lh[l + 1];
    auto n = end - begin;
    selection_sort(SortableAdj{lh2h, codes, hg, begin}, n);
  //CHECK(n <= MAX_UPWARD);
  //top_down_merge_sort<MAX_UPWARD>(SortableAdj{lh2h, codes, hg, begin}, n);
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
