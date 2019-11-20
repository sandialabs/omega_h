#include "Omega_h_swap3d.hpp"

#include "Omega_h_for.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_swap3d_choice.hpp"
#include "Omega_h_swap3d_loop.hpp"

namespace Omega_h {

template <Int metric_dim>
void swap3d_qualities_tmpl(Mesh* mesh, AdaptOpts const& opts,
    LOs cands2edges, Reals* cand_quals, Read<I8>* cand_configs) {
  auto edges2tets = mesh->ask_up(EDGE, REGION);
  auto edges2edge_tets = edges2tets.a2ab;
  auto edge_tets2tets = edges2tets.ab2b;
  auto edge_tet_codes = edges2tets.codes;
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto tet_verts2verts = mesh->ask_verts_of(REGION);
  auto edges_are_owned = mesh->owned(EDGE);
  auto quality_measure = MetricElementQualities<3, metric_dim>(mesh);
  auto length_measure = MetricEdgeLengths<3, metric_dim>(mesh);
  auto max_length = opts.max_length_allowed;
  auto ncands = cands2edges.size();
  auto cand_quals_w = Write<Real>(ncands);
  auto cand_configs_w = Write<I8>(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto edge = cands2edges[cand];
    /* non-owned edges will have incomplete cavities
       and will run into the topological assertions
       in find_loop(). don't bother; their results
       will be overwritten by the owner's anyways */
    if (!edges_are_owned[edge]) {
      cand_configs_w[cand] = -1;
      cand_quals_w[cand] = -1.0;
      return;
    }
    auto loop = swap3d::find_loop(edges2edge_tets, edge_tets2tets,
        edge_tet_codes, edge_verts2verts, tet_verts2verts, edge);
    if (loop.size > swap3d::MAX_EDGE_SWAP) {
      cand_configs_w[cand] = -1;
      cand_quals_w[cand] = -1.0;
      return;
    }
    auto choice =
        swap3d::choose(loop, quality_measure, length_measure, max_length);
    static_assert(swap3d::MAX_CONFIGS <= INT8_MAX,
        "int8_t must be able to represent all swap configurations");
    cand_configs_w[cand] = static_cast<I8>(choice.mesh);
    cand_quals_w[cand] = choice.quality;
  };
  parallel_for(ncands, f, "swap3d_qualities");
  *cand_quals = cand_quals_w;
  *cand_configs = cand_configs_w;
  *cand_quals =
      mesh->sync_subset_array(EDGE, *cand_quals, cands2edges, -1.0, 1);
  *cand_configs =
      mesh->sync_subset_array(EDGE, *cand_configs, cands2edges, I8(-1), 1);
}

void swap3d_qualities(Mesh* mesh, AdaptOpts const& opts, LOs cands2edges,
    Reals* cand_quals, Read<I8>* cand_configs) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  OMEGA_H_CHECK(mesh->dim() == 3);
  auto metrics = mesh->get_array<Real>(VERT, "metric");
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (metric_dim == 3) {
    swap3d_qualities_tmpl<3>(mesh, opts, cands2edges, cand_quals, cand_configs);
    return;
  }
  if (metric_dim == 1) {
    swap3d_qualities_tmpl<1>(mesh, opts, cands2edges, cand_quals, cand_configs);
    return;
  }
  OMEGA_H_NORETURN();
}

}  // end namespace Omega_h
