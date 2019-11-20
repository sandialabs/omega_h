#include "Omega_h_swap2d.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

/* for similarity to swap3d and to take advantage
   of the existing computation of the two relevant
   vertices, edge length overshooting is also handled
   by the qualities function */

template <Int metric_dim>
Reals swap2d_qualities_tmpl(
    Mesh* mesh, AdaptOpts const& opts, LOs cands2edges) {
  auto e2t = mesh->ask_up(EDGE, FACE);
  auto e2et = e2t.a2ab;
  auto et2t = e2t.ab2b;
  auto et_codes = e2t.codes;
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto tri_verts2verts = mesh->ask_verts_of(FACE);
  auto edges_are_owned = mesh->owned(EDGE);
  auto quality_measure = MetricElementQualities<2, metric_dim>(mesh);
  auto length_measure = MetricEdgeLengths<2, metric_dim>(mesh);
  auto max_length = opts.max_length_allowed;
  auto ncands = cands2edges.size();
  auto cand_quals_w = Write<Real>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto edge = cands2edges[cand];
    /* non-owned edges will have incomplete cavities
       and will run into the topological assertions
       in find_loop(). don't bother; their results
       will be overwritten by the owner's anyways */
    if (!edges_are_owned[edge]) {
      cand_quals_w[cand] = -1.0;
      return;
    }
    OMEGA_H_CHECK(e2et[edge + 1] == 2 + e2et[edge]);
    LO t[2];
    Few<LO, 2> ov;
    for (Int i = 0; i < 2; ++i) {
      auto et = e2et[edge] + i;
      auto code = et_codes[et];
      auto tte = code_which_down(code);
      auto rot = code_rotation(code);
      t[rot] = et2t[et];
      auto ttv = simplex_opposite_template(FACE, EDGE, tte);
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
  parallel_for(ncands, f, "swap2d_qualities");
  auto cand_quals = Reals(cand_quals_w);
  return mesh->sync_subset_array(EDGE, cand_quals, cands2edges, -1.0, 1);
}

Reals swap2d_qualities(Mesh* mesh, AdaptOpts const& opts, LOs cands2edges) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto metrics = mesh->get_array<Real>(VERT, "metric");
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (metric_dim == 2) {
    return swap2d_qualities_tmpl<2>(mesh, opts, cands2edges);
  }
  if (metric_dim == 1) {
    return swap2d_qualities_tmpl<1>(mesh, opts, cands2edges);
  }
  OMEGA_H_NORETURN(Reals());
}

}  // end namespace Omega_h
