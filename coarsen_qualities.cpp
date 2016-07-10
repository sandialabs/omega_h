#include "internal.hpp"

namespace osh {

template <typename Measure, Int dim>
static Reals coarsen_qualities_tmpl(Mesh* mesh, LOs cands2edges,
                                    Read<I8> cand_codes) {
  Measure measure(mesh);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto cv2v = mesh->ask_verts_of(dim);
  auto v2c = mesh->ask_up(VERT, dim);
  auto v2vc = v2c.a2ab;
  auto vc2c = v2c.ab2b;
  auto vc_codes = v2c.codes;
  auto ncands = cands2edges.size();
  auto qualities = Write<Real>(ncands * 2, -1.0);
  auto f = LAMBDA(LO cand) {
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
        auto ccv2v = gather_verts<dim + 1>(cv2v, c);
        bool will_die = false;
        for (auto ccv = 0; ccv < (dim + 1); ++ccv) {
          if ((ccv != ccv_col) && (ccv2v[ccv] == v_onto)) {
            will_die = true;
            break;
          }
        }
        if (will_die) continue;
        ccv2v[ccv_col] = v_onto;  // vertices of new cell
        auto qual = measure.measure(ccv2v);
        minqual = min2(minqual, qual);
      }
      qualities[cand * 2 + eev_col] = minqual;
    }
  };
  parallel_for(ncands, f);
  return qualities;
}

Reals coarsen_qualities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes) {
  CHECK(mesh->parting() == OSH_GHOSTED);
  auto cand_quals = Reals();
  if (mesh->dim() == 3) {
    if (mesh->has_tag(VERT, "metric")) {
      cand_quals = coarsen_qualities_tmpl<MetricElementQualities, 3>(
          mesh, cands2edges, cand_codes);
    } else {
      cand_quals = coarsen_qualities_tmpl<RealElementQualities, 3>(
          mesh, cands2edges, cand_codes);
    }
  } else {
    CHECK(mesh->dim() == 2);
    if (mesh->has_tag(VERT, "metric")) {
      cand_quals = coarsen_qualities_tmpl<MetricElementQualities, 2>(
          mesh, cands2edges, cand_codes);
    } else {
      cand_quals = coarsen_qualities_tmpl<RealElementQualities, 2>(
          mesh, cands2edges, cand_codes);
    }
  }
  return mesh->sync_subset_array(EDGE, cand_quals, cands2edges, -1.0, 2);
}

void choose_vertex_collapses(Mesh* mesh, LOs cands2edges,
                             Read<I8> cand_edge_codes, Reals cand_edge_quals,
                             Read<I8>& verts_are_cands, Reals& vert_quals) {
  CHECK(mesh->parting() == OSH_GHOSTED);
  auto edges2cands = invert_injective_map(cands2edges, mesh->nedges());
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto ve_codes = v2e.codes;
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto verts_are_cands_w = Write<I8>(mesh->nverts());
  auto vert_quals_w = Write<Real>(mesh->nverts());
  auto f = LAMBDA(LO v) {
    bool vert_is_cand = false;
    Real best_qual = -1.0;
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto cand = edges2cands[e];
      if (cand == -1) continue;
      auto ve_code = ve_codes[ve];
      auto eev = code_which_down(ve_code);
      auto cand_code = cand_edge_codes[cand];
      if (!collapses(cand_code, eev)) continue;
      auto qual = cand_edge_quals[cand * 2 + eev];
      if (qual > best_qual) {
        vert_is_cand = true;
        best_qual = qual;
      }
    }
    verts_are_cands_w[v] = vert_is_cand;
    vert_quals_w[v] = best_qual;
  };
  parallel_for(mesh->nverts(), f);
  verts_are_cands = verts_are_cands_w;
  verts_are_cands = mesh->sync_array(VERT, verts_are_cands, 1);
  vert_quals = vert_quals_w;
  vert_quals = mesh->sync_array(VERT, vert_quals, 1);
}

static Read<I8> filter_coarsen_dirs(Read<I8> codes, Read<I8> keep_dirs) {
  auto codes_w = Write<I8>(codes.size());
  auto f = LAMBDA(LO cand) {
    auto code = codes[cand];
    for (Int dir = 0; dir < 2; ++dir) {
      if (!keep_dirs[cand * 2 + dir]) {
        code = dont_collapse(code, dir);
      }
    }
    codes_w[cand] = code;
  };
  parallel_for(codes_w.size(), f);
  return codes_w;
}

Read<I8> filter_coarsen_min_qual(Read<I8> cand_codes, Reals cand_quals,
                                 Real min_qual) {
  auto keep_dirs = each_geq_to(cand_quals, min_qual);
  return filter_coarsen_dirs(cand_codes, keep_dirs);
}

Read<I8> filter_coarsen_improve(Mesh* mesh, LOs cands2edges,
                                Read<I8> cand_codes, Reals cand_quals) {
  auto elem_quals = mesh->ask_qualities();
  auto verts2elems = mesh->ask_up(VERT, mesh->dim());
  auto vert_old_quals = graph_reduce(verts2elems, elem_quals, 1, OSH_MIN);
  vert_old_quals = mesh->sync_array(VERT, vert_old_quals, 1);
  auto edge_verts2verts = mesh->ask_verts_of(EDGE);
  auto edge_old_quals = unmap(edge_verts2verts, vert_old_quals, 1);
  auto cand_old_quals = unmap(cands2edges, edge_old_quals, 2);
  auto keep_dirs = gt_each(cand_quals, cand_old_quals);
  return filter_coarsen_dirs(cand_codes, keep_dirs);
}

} //end namespace osh
