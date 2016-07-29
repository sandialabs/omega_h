#include "adjacency.hpp"

#include "algorithm.hpp"
#include "align.hpp"
#include "local_sort.hpp"
#include "loop.hpp"
#include "map.hpp"

#define MAX_UPWARD 100

namespace osh {

static void order_by_globals(LOs l2lh, Write<LO> lh2h, Write<I8> codes,
                             Read<GO> hg) {
  struct Item {
    GO high_global;
    LO high;
    I8 code;
    DEVICE bool operator<(Item const& other) {
      return high_global < other.high_global;
    }
  };
  LO nl = l2lh.size() - 1;
  auto f = LAMBDA(LO l) {
    LO begin = l2lh[l];
    LO end = l2lh[l + 1];
    CHECK(end - begin <= MAX_UPWARD);
    LocalMergeSort<Item, MAX_UPWARD> sorter;
    sorter.n = end - begin;
    for (Int i = 0; i < sorter.n; ++i) {
      auto lh = begin + i;
      Item item;
      item.high = lh2h[lh];
      item.code = codes[lh];
      item.high_global = hg[item.high];
      sorter.array[i] = item;
    }
    sorter.run();
    for (Int i = 0; i < sorter.n; ++i) {
      auto lh = begin + i;
      auto item = sorter.array[i];
      lh2h[lh] = item.high;
      codes[lh] = item.code;
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
