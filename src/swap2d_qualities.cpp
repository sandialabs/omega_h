#include "swap2d.hpp"

#include "access.hpp"
#include "quality.hpp"
#include "simplices.hpp"

namespace Omega_h {

/* for similarity to swap3d and to take advantage
   of the existing computation of the two relevant
   vertices, edge length overshooting is also handled
   by the qualities function */

template <typename QualityMeasure, typename LengthMeasure>
static Reals swap2d_qualities_tmpl(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2edges) {
  auto e2t = mesh->ask_up(EDGE, TRI);
  auto e2et = e2t.a2ab;
  auto et2t = e2t.ab2b;
  auto et_codes = e2t.codes;
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto tri_verts2verts = mesh->ask_verts_of(TRI);
  auto edges_are_owned = mesh->owned(EDGE);
  auto quality_measure = QualityMeasure(mesh);
  auto length_measure = LengthMeasure(mesh);
  auto max_length = opts.max_length_allowed;
  auto ncands = cands2edges.size();
  auto cand_quals_w = Write<Real>(ncands);
  auto f = LAMBDA(LO cand) {
    auto edge = cands2edges[cand];
    /* non-owned edges will have incomplete cavities
       and will run into the topological assertions
       in find_loop(). don't bother; their results
       will be overwritten by the owner's anyways */
    if (!edges_are_owned[edge]) {
      cand_quals_w[cand] = -1.0;
      return;
    }
    CHECK(e2et[edge + 1] == 2 + e2et[edge]);
    LO t[2];
    Few<LO, 2> ov;
    for (Int i = 0; i < 2; ++i) {
      auto et = e2et[edge] + i;
      auto code = et_codes[et];
      auto tte = code_which_down(code);
      auto rot = code_rotation(code);
      t[rot] = et2t[et];
      auto ttv = opposite_template(TRI, EDGE, tte);
      ov[rot] = tri_verts2verts[t[rot] * 3 + ttv];
    }
    auto l = length_measure.measure(ov);
    if (l > max_length) {
      cand_quals_w[cand] = -1.0;
      return;
    }
    auto ev = gather_verts<2>(edge_verts2verts, edge);
    Real minqual = 1.0;
    for (Int i = 0; i < 2; ++i) {
      Few<LO, 3> ntv;
      ntv[0] = ev[1 - i];
      ntv[1] = ov[i];
      ntv[2] = ov[1 - i];
      auto qual = quality_measure.measure(ntv);
      minqual = min2(minqual, qual);
    }
    cand_quals_w[cand] = minqual;
  };
  parallel_for(ncands, f);
  return cand_quals_w;
}

Reals swap2d_qualities(Mesh* mesh, AdaptOpts const& opts, LOs cands2edges) {
  CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto cand_quals = Reals();
  if (mesh->has_tag(VERT, "metric")) {
    using MEQ = MetricElementQualities;
    using MEL = MetricEdgeLengths<2>;
    cand_quals = swap2d_qualities_tmpl<MEQ, MEL>(mesh, opts, cands2edges);
  } else {
    CHECK(mesh->has_tag(VERT, "size"));
    using REQ = RealElementQualities;
    using REL = IsoEdgeLengths<2>;
    cand_quals = swap2d_qualities_tmpl<REQ, REL>(mesh, opts, cands2edges);
  }
  return mesh->sync_subset_array(EDGE, cand_quals, cands2edges, -1.0, 1);
}

}  // end namespace Omega_h
