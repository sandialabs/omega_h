#include "Omega_h_swap2d.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

void swap2d_topology(Mesh* mesh, LOs keys2edges,
    HostFew<LOs, 3>* keys2prods_out, HostFew<LOs, 3>* prod_verts2verts_out) {
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto tv2v = mesh->ask_verts_of(FACE);
  auto e2t = mesh->ask_up(EDGE, FACE);
  auto e2et = e2t.a2ab;
  auto et2t = e2t.ab2b;
  auto et_codes = e2t.codes;
  auto nkeys = keys2edges.size();
  auto tri_verts2verts_w = Write<LO>(nkeys * 2 * 3);
  auto edge_verts2verts_w = Write<LO>(nkeys * 1 * 2);
  auto f = OMEGA_H_LAMBDA(LO key) {
    auto e = keys2edges[key];
    OMEGA_H_CHECK(e2et[e + 1] == 2 + e2et[e]);
    LO t[2];
    auto ov = &edge_verts2verts_w[key * 2];
    for (Int i = 0; i < 2; ++i) {
      auto et = e2et[e] + i;
      auto code = et_codes[et];
      auto tte = code_which_down(code);
      auto rot = code_rotation(code);
      t[rot] = et2t[et];
      auto ttv = simplex_opposite_template(FACE, EDGE, tte);
      ov[rot] = tv2v[t[rot] * 3 + ttv];
    }
    auto ev = gather_verts<2>(ev2v, e);
    for (Int i = 0; i < 2; ++i) {
      auto prod = key * 2 + i;
      auto ntv = &tri_verts2verts_w[prod * 3];
      ntv[0] = ev[1 - i];
      ntv[1] = ov[i];
      ntv[2] = ov[1 - i];
    }
  };
  parallel_for(nkeys, f, "swap2d_topology");
  HostFew<LOs, 3> keys2prods;
  HostFew<LOs, 3> prod_verts2verts;
  keys2prods[EDGE] = LOs(nkeys + 1, 0, 1);
  keys2prods[FACE] = LOs(nkeys + 1, 0, 2);
  prod_verts2verts[EDGE] = edge_verts2verts_w;
  prod_verts2verts[FACE] = tri_verts2verts_w;
  *keys2prods_out = keys2prods;
  *prod_verts2verts_out = prod_verts2verts;
}

}  // end namespace Omega_h
