template <typename Measure, Int dim>
static Reals coarsen_qualities_tmpl(Mesh& mesh,
    LOs cands2edges,
    Read<I8> cand_codes) {
  Measure measure(mesh);
  auto ev2v = mesh.ask_verts_of(EDGE);
  auto cv2v = mesh.ask_verts_of(dim);
  auto v2c = mesh.ask_up(VERT, dim);
  auto v2vc = v2c.a2ab;
  auto vc2c = v2c.ab2b;
  auto vc_codes = v2c.codes;
  auto ncands = cands2edges.size();
  auto qualities = Write<Real>(ncands * 2);
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
        ccv2v[ccv_col] = v_onto; //vertices of new cell
        auto qual = measure.measure(ccv2v);
        minqual = min2(minqual, qual);
      }
      qualities[cand * 2 + eev_col] = minqual;
    }
  };
  parallel_for(ncands, f);
  return qualities;
}

Reals coarsen_qualities(Mesh& mesh,
    LOs cands2edges,
    Read<I8> cand_codes) {
  CHECK(mesh.partition() == GHOSTED);
  if (mesh.dim() == 3) {
    if (mesh.has_tag(VERT, "metric")) {
      return coarsen_qualities_tmpl<MetricElementQualities,3>(
          mesh, cands2edges, cand_codes);
    } else {
      return coarsen_qualities_tmpl<RealElementQualities,3>(
          mesh, cands2edges, cand_codes);
    }
  } else {
    CHECK(mesh.dim() == 2);
    if (mesh.has_tag(VERT, "metric")) {
      return coarsen_qualities_tmpl<MetricElementQualities,2>(
          mesh, cands2edges, cand_codes);
    } else {
      return coarsen_qualities_tmpl<RealElementQualities,2>(
          mesh, cands2edges, cand_codes);
    }
  }
}
