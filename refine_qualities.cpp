struct RealRefineQualities {
  RealRefineQualities(Mesh const&) {}
  template <Int dim>
  struct Helper {
  };
  template <Int dim>
  Helper<dim> get_helper(Few<LO, 2>) const {
    return Helper<dim>();
  }
  template <Int dim>
  INLINE Real measure(Few<Vector<dim>, dim + 1> p,
      Helper<dim>) const {
    return real_element_quality(p);
  }
};

struct MetricRefineQualities {
  Reals metrics;
  MetricRefineQualities(Mesh const& mesh):
    metrics(mesh.get_array<Real>(VERT, "metric"))
  {}
  template <Int dim>
  struct Helper {
    Matrix<dim,dim> metric;
  };
  template <Int dim>
  Helper<dim> get_helper(Few<LO, 2> eev2v) const {
    Helper<dim> helper;
    helper.metric = average_metrics(gather_symms<2,dim>(metrics, eev2v));
    return helper;
  }
  template <Int dim>
  INLINE Real measure(Few<Vector<dim>, dim + 1> p,
      Helper<dim> helper) const {
    return metric_element_quality(p, helper.metric);
  }
};

template <typename Measure, Int dim>
static Reals refine_qualities_tmpl(Mesh& mesh, LOs candidates) {
  auto ev2v = mesh.ask_down(EDGE, VERT).ab2b;
  auto cv2v = mesh.ask_down(dim, VERT).ab2b;
  auto e2c = mesh.ask_up(EDGE, dim);
  auto e2ec = e2c.a2ab;
  auto ec2c = e2c.ab2b;
  auto ec_codes = e2c.codes;
  auto coords = mesh.coords();
  auto ncands = candidates.size();
  auto measure = Measure(mesh);
  Write<Real> quals(ncands);
  auto f = LAMBDA(LO cand) {
    auto e = candidates[cand];
    auto eev2v = gather_verts<2>(ev2v, e);
    typename Measure::template Helper<dim> helper;
    helper = measure.template get_helper<dim>(eev2v);
    auto ep = gather_vectors<2, dim>(coords, eev2v);
    auto midp = (ep[0] + ep[1]) / 2.;
    auto minqual = 1.0;
    for (auto ec = e2ec[e]; ec < e2ec[e + 1]; ++ec) {
      auto c = ec2c[ec];
      auto code = ec_codes[ec];
      auto cce = code_which_down(code);
      auto rot = code_rotation(code);
      auto ccv2v = gather_verts<dim + 1>(cv2v, c);
      for (Int eev = 0; eev < 2; ++eev) {
        /* a new cell is formed from an old cell by finding
           its side that is opposite to one of the edge endpoints
           and connecting it to the midpoint to form the new cell */
        auto cev = eev ^ rot;
        auto ccv = down_templates[dim][EDGE][cce][cev];
        auto ccs = opposite_templates[dim][VERT][ccv];
        Few<Vector<dim>, dim + 1> ncp;
        for (Int csv = 0; csv < dim; ++csv) {
          auto ccv2 = down_templates[dim][dim - 1][ccs][csv];
          auto v2 = ccv2v[ccv2];
          ncp[csv] = get_vec<dim>(coords, v2);
        }
        ncp[dim] = midp;
        auto cqual = measure.measure(ncp, helper);
        minqual = min2(minqual, cqual);
      }
    }
    quals[cand] = minqual;
  };
  parallel_for(ncands, f);
  return quals;
}

Reals refine_qualities(Mesh& mesh, LOs candidates) {
  auto dim = mesh.dim();
  auto have_metric = mesh.has_tag(VERT, "metric");
  if (have_metric) {
    if (dim == 3) {
      return refine_qualities_tmpl<MetricRefineQualities,3>(
          mesh, candidates);
    } else {
      CHECK(dim == 2);
      return refine_qualities_tmpl<MetricRefineQualities,2>(
          mesh, candidates);
    }
  } else {
    if (dim == 3) {
      return refine_qualities_tmpl<RealRefineQualities,3>(
          mesh, candidates);
    } else {
      CHECK(dim == 2);
      return refine_qualities_tmpl<RealRefineQualities,2>(
          mesh, candidates);
    }
  }
}

