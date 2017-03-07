#include "adjacency.hpp"

#include "internal.hpp"
#include "map.hpp"
#include "scan.hpp"
#include "loop.hpp"

namespace Omega_h {

Adj unmap_adjacency(LOs a2b, Adj b2c) {
  auto b2bc = b2c.a2ab;
  auto bc2c = b2c.ab2b;
  auto bc_codes = b2c.codes;
  auto b_degrees = get_degrees(b2bc);
  auto a_degrees = unmap(a2b, b_degrees, 1);
  auto a2ac = offset_scan(a_degrees);
  auto na = a2b.size();
  auto nac = a2ac.last();
  Write<LO> ac2c(nac);
  auto ac_codes = Write<I8>(nac);
  auto f = LAMBDA(LO a) {
    auto b = a2b[a];
    auto bc = b2bc[b];
    for (auto ac = a2ac[a]; ac < a2ac[a + 1]; ++ac) {
      ac2c[ac] = bc2c[bc];
      ac_codes[ac] = bc_codes[bc];
      ++bc;
    }
  };
  parallel_for(na, f);
  return Adj(a2ac, ac2c, ac_codes);
}

}  // namespace Omega_h
