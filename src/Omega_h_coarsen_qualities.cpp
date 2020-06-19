#include "Omega_h_coarsen.hpp"

#include "Omega_h_align.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_collapse.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_quality.hpp"

#include <iostream>

namespace Omega_h {

template <Int mesh_dim, Int metric_dim>
Reals coarsen_qualities_tmpl(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  OMEGA_H_CHECK(mesh->dim() == mesh_dim);
  MetricElementQualities<mesh_dim, metric_dim> measure(mesh);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto cv2v = mesh->ask_elem_verts();
  auto v2c = mesh->ask_up(VERT, mesh_dim);
  auto v2vc = v2c.a2ab;
  auto vc2c = v2c.ab2b;
  auto vc_codes = v2c.codes;
  auto ncands = cands2edges.size();
  auto qualities = Write<Real>(ncands * 2, -1.0);
  auto coords = mesh->coords();
  auto is_bad_w = Write<Byte>(mesh->nelems(), Byte(0));
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto e = cands2edges[cand];
    auto code = cand_codes[cand];
    for (Int eev_col = 0; eev_col < 2; ++eev_col) {
      if (!collapses(code, eev_col)) continue;
      auto v_col = ev2v[e * 2 + eev_col];
      auto eev_onto = 1 - eev_col;
      auto v_onto = ev2v[e * 2 + eev_onto];
      Real minqual = 1.0;
      for (auto vc = v2vc[v_col]; vc < v2vc[v_col + 1]; ++vc) {
        auto c = vc2c[vc];
        auto vc_code = vc_codes[vc];
        auto ccv_col = code_which_down(vc_code);
        auto ccv2v = gather_verts<mesh_dim + 1>(cv2v, c);
        bool will_die = false;
        for (auto ccv = 0; ccv < (mesh_dim + 1); ++ccv) {
          if ((ccv != ccv_col) && (ccv2v[ccv] == v_onto)) {
            will_die = true;
            break;
          }
        }
        if (will_die) continue;
        OMEGA_H_CHECK(0 <= ccv_col && ccv_col < mesh_dim + 1);
        ccv2v[ccv_col] = v_onto;  // vertices of new cell
        auto qual = measure.measure(ccv2v);
        minqual = min2(minqual, qual);
      }
      qualities[cand * 2 + eev_col] = minqual;
    }
  };
  parallel_for(ncands, f, "coarsen_qualities");
  auto out = Reals(qualities);
  return mesh->sync_subset_array(EDGE, out, cands2edges, -1.0, 2);
}

Reals coarsen_qualities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  OMEGA_H_CHECK(mesh->parting() == OMEGA_H_GHOSTED);
  auto metrics = mesh->get_array<Real>(VERT, "metric");
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  auto cand_quals = Reals();
  if (mesh->dim() == 3 && metric_dim == 3) {
    return coarsen_qualities_tmpl<3, 3>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return coarsen_qualities_tmpl<2, 2>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return coarsen_qualities_tmpl<3, 1>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return coarsen_qualities_tmpl<2, 1>(mesh, cands2edges, cand_codes);
  }
  if (mesh->dim() == 1) {
    auto edges2verts = mesh->ask_verts_of(EDGE);
    auto cands2verts = unmap(cands2edges, edges2verts, 2);
    return get_1d_cavity_qualities(mesh, VERT, cands2verts);
  }
  OMEGA_H_NORETURN(Reals());
}

Read<I8> filter_coarsen_dirs(Read<I8> codes, Read<I8> keep_dirs) {
  auto codes_w = Write<I8>(codes.size());
  auto f = OMEGA_H_LAMBDA(LO cand) {
    auto code = codes[cand];
    for (Int dir = 0; dir < 2; ++dir) {
      if (!keep_dirs[cand * 2 + dir]) {
        code = dont_collapse(code, dir);
      }
    }
    codes_w[cand] = code;
  };
  parallel_for(codes_w.size(), f, "filter_coarsen_dirs");
  return codes_w;
}

Read<I8> filter_coarsen_min_qual(
    Read<I8> cand_codes, Reals cand_quals, Real min_qual) {
  auto keep_dirs = each_geq_to(cand_quals, min_qual);
  return filter_coarsen_dirs(cand_codes, keep_dirs);
}

Read<I8> filter_coarsen_improve(
    Mesh* mesh, LOs cands2edges, Read<I8> cand_codes, Reals cand_quals) {
  auto elem_quals = mesh->ask_qualities();
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto vert_old_quals = graph_reduce(verts2elems, elem_quals, 1, OMEGA_H_MIN);
  vert_old_quals = mesh->sync_array(VERT, vert_old_quals, 1);
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto edge_old_quals = read(unmap(edge_verts2verts, vert_old_quals, 1));
  auto cand_old_quals = read(unmap(cands2edges, edge_old_quals, 2));
  auto keep_dirs = gt_each(cand_quals, cand_old_quals);
  return filter_coarsen_dirs(cand_codes, keep_dirs);
}

}  // end namespace Omega_h
