#include "swap2d.hpp"

#include "access.hpp"
#include "loop.hpp"
#include "quality.hpp"
#include "simplices.hpp"

namespace Omega_h {

template <typename Measure, Int dim>
static Reals swap2d_qualities_tmpl(Mesh* mesh, LOs cands2edges) {
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto tv2v = mesh->ask_verts_of(TRI);
  auto e2t = mesh->ask_up(EDGE, TRI);
  auto e2et = e2t.a2ab;
  auto et2t = e2t.ab2b;
  auto et_codes = e2t.codes;
  auto measure = Measure(mesh);
  auto ncands = cands2edges.size();
  auto quals_w = Write<Real>(ncands);
  auto f = LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    CHECK(e2et[e + 1] == 2 + e2et[e]);
    LO t[2];
    LO ov[2];
    for (Int i = 0; i < 2; ++i) {
      auto et = e2et[e] + i;
      auto code = et_codes[et];
      auto tte = code_which_down(code);
      auto rot = code_rotation(code);
      t[rot] = et2t[et];
      auto ttv = OppositeTemplate<TRI, EDGE>::get(tte);
      ov[rot] = tv2v[t[rot] * 3 + ttv];
    }
    auto ev = gather_verts<2>(ev2v, e);
    Real minqual = 1.0;
    for (Int i = 0; i < 2; ++i) {
      Few<LO, 3> ntv;
      ntv[0] = ev[1 - i];
      ntv[1] = ov[i];
      ntv[2] = ov[1 - i];
      auto qual = measure.measure(ntv);
      minqual = min2(minqual, qual);
    }
    quals_w[cand] = minqual;
  };
  parallel_for(ncands, f);
  return quals_w;
}

Reals swap2d_qualities(Mesh* mesh, LOs cands2edges) {
  CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto cand_quals = Reals();
  if (mesh->dim() == 3) {
    if (mesh->has_tag(VERT, "metric")) {
      cand_quals =
          swap2d_qualities_tmpl<MetricElementQualities, 3>(mesh, cands2edges);
    } else {
      cand_quals =
          swap2d_qualities_tmpl<RealElementQualities, 3>(mesh, cands2edges);
    }
  } else {
    CHECK(mesh->dim() == 2);
    if (mesh->has_tag(VERT, "metric")) {
      cand_quals =
          swap2d_qualities_tmpl<MetricElementQualities, 2>(mesh, cands2edges);
    } else {
      cand_quals =
          swap2d_qualities_tmpl<RealElementQualities, 2>(mesh, cands2edges);
    }
  }
  return mesh->sync_subset_array(EDGE, cand_quals, cands2edges, -1.0, 1);
}

}  // end namespace Omega_h
