#include "Omega_h_refine_qualities.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_quality.hpp"
#include "Omega_h_refine_topology.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
struct MetricRefineQualities {
  Reals vert_metrics;
  Reals midpt_metrics;
  MetricRefineQualities(Mesh* mesh, LOs candidates)
      : vert_metrics(mesh->get_array<Real>(VERT, "metric")),
        /* TODO: we could reuse the results of this instead of recomputing
         * them when transferring an OMEGA_H_METRIC field in transfer.cpp
         */
        midpt_metrics(get_mident_metrics(
            mesh, EDGE, candidates, mesh->get_array<Real>(VERT, "metric"))) {}
  OMEGA_H_DEVICE Real measure(Int cand, Few<Vector<mesh_dim>, mesh_dim + 1> p,
      Few<LO, mesh_dim> csv2v) const {
    Few<Matrix<metric_dim, metric_dim>, mesh_dim + 1> ms;
    for (Int csv = 0; csv < mesh_dim; ++csv)
      ms[csv] = get_symm<metric_dim>(vert_metrics, csv2v[csv]);
    ms[mesh_dim] = get_symm<metric_dim>(midpt_metrics, cand);
    auto m = maxdet_metric(ms);
    return metric_element_quality(p, m);
  }
};

template <Int mesh_dim, Int metric_dim>
Reals refine_qualities_tmpl(Mesh* mesh, LOs candidates) {
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto cv2v = mesh->ask_verts_of(mesh_dim);
  auto e2c = mesh->ask_up(EDGE, mesh_dim);
  auto e2ec = e2c.a2ab;
  auto ec2c = e2c.ab2b;
  auto ec_codes = e2c.codes;
  auto coords = mesh->coords();
  auto ncands = candidates.size();
  auto measure = MetricRefineQualities<mesh_dim, metric_dim>(mesh, candidates);
  Write<Real> quals_w(ncands);
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto e = candidates[cand];
    auto eev2v = gather_verts<2>(ev2v, e);
    auto ep = gather_vectors<2, mesh_dim>(coords, eev2v);
    auto midp = (ep[0] + ep[1]) / 2.;
    auto minqual = 1.0;
    for (auto ec = e2ec[e]; ec < e2ec[e + 1]; ++ec) {
      auto c = ec2c[ec];
      auto code = ec_codes[ec];
      auto cce = code_which_down(code);
      auto rot = code_rotation(code);
      auto ccv2v = gather_verts<mesh_dim + 1>(cv2v, c);
      for (Int eev = 0; eev < 2; ++eev) {
        /* a new cell is formed from an old cell by finding
           its side that is opposite to one of the edge endpoints
           and connecting it to the midpoint to form the new cell
           (see refine_domain_interiors) */
        auto cev = eev ^ rot;
        auto ccv = simplex_down_template(mesh_dim, EDGE, cce, cev);
        auto ccs = simplex_opposite_template(mesh_dim, VERT, ccv);
        Few<LO, mesh_dim> csv2v;
        Few<Vector<mesh_dim>, mesh_dim + 1> ncp;
        for (Int csv = 0; csv < mesh_dim; ++csv) {
          auto ccv2 = simplex_down_template(mesh_dim, mesh_dim - 1, ccs, csv);
          auto v2 = ccv2v[ccv2];
          csv2v[csv] = v2;
          ncp[csv] = get_vector<mesh_dim>(coords, v2);
        }
        ncp[mesh_dim] = midp;
        flip_new_elem<mesh_dim>(&csv2v[0]);
        flip_new_elem<mesh_dim>(&ncp[0]);
        auto cqual = measure.measure(cand, ncp, csv2v);
        minqual = min2(minqual, cqual);
      }
    }
    quals_w[cand] = minqual;
  };
  parallel_for(ncands, f, "refine_qualities");
  auto cand_quals = Reals(quals_w);
  return mesh->sync_subset_array(EDGE, cand_quals, candidates, -1.0, 1);
}

Reals refine_qualities(Mesh* mesh, LOs candidates) {
  auto mesh_dim = mesh->dim();
  auto metric_dim = get_metric_dim(mesh);
  if (mesh_dim == 3 && metric_dim == 3) {
    return refine_qualities_tmpl<3, 3>(mesh, candidates);
  }
  if (mesh_dim == 2 && metric_dim == 2) {
    return refine_qualities_tmpl<2, 2>(mesh, candidates);
  }
  if (mesh_dim == 3 && metric_dim == 1) {
    return refine_qualities_tmpl<3, 1>(mesh, candidates);
  }
  if (mesh_dim == 2 && metric_dim == 1) {
    return refine_qualities_tmpl<2, 1>(mesh, candidates);
  }
  if (mesh_dim == 1) {
    return get_1d_cavity_qualities(mesh, EDGE, candidates);
  }
  OMEGA_H_NORETURN(Reals());
}

}  // end namespace Omega_h
