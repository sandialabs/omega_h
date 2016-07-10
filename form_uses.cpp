#include "adjacency.hpp"

#include "loop.hpp"
#include "simplices.hpp"

namespace osh {

LOs form_uses(LOs hv2v, Int high_dim, Int low_dim) {
  Int nverts_per_high = simplex_degrees[high_dim][0];
  Int nverts_per_low = simplex_degrees[low_dim][0];
  Int nlows_per_high = simplex_degrees[high_dim][low_dim];
  LO nhigh = hv2v.size() / nverts_per_high;
  LO nuses = nhigh * nlows_per_high;
  Write<LO> uv2v(nuses * nverts_per_low);
  auto f = LAMBDA(LO h) {
    LO h_begin = h * nverts_per_high;
    for (Int u = 0; u < nlows_per_high; ++u) {
      LO u_begin = (h * nlows_per_high + u) * nverts_per_low;
      for (Int uv = 0; uv < nverts_per_low; ++uv) {
        uv2v[u_begin + uv] =
            hv2v[h_begin + down_templates[high_dim][low_dim][u][uv]];
      }
    }
  };
  parallel_for(nhigh, f);
  return uv2v;
}

} //end namespace osh
