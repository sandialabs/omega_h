#include "size.hpp"

#include "array.hpp"
#include "graph.hpp"
#include "loop.hpp"

namespace osh {

template <typename EdgeLengths>
static Reals measure_edges_tmpl(Mesh* mesh, LOs a2e) {
  EdgeLengths measurer(mesh);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto na = a2e.size();
  Write<Real> lengths(na);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<2>(ev2v, e);
    lengths[a] = measurer.measure(v);
  };
  parallel_for(na, f);
  return lengths;
}

Reals measure_edges_real(Mesh* mesh, LOs a2e) {
  if (mesh->dim() == 3) {
    return measure_edges_tmpl<RealEdgeLengths<3>>(mesh, a2e);
  } else {
    CHECK(mesh->dim() == 2);
    return measure_edges_tmpl<RealEdgeLengths<2>>(mesh, a2e);
  }
}

Reals measure_edges_metric(Mesh* mesh, LOs a2e) {
  if (mesh->dim() == 3) {
    if (mesh->has_tag(VERT, "size")) {
      return measure_edges_tmpl<IsoEdgeLengths<3>>(mesh, a2e);
    }
    if (mesh->has_tag(VERT, "metric")) {
      return measure_edges_tmpl<MetricEdgeLengths<3>>(mesh, a2e);
    }
  } else {
    CHECK(mesh->dim() == 2);
    if (mesh->has_tag(VERT, "size")) {
      return measure_edges_tmpl<IsoEdgeLengths<2>>(mesh, a2e);
    }
    if (mesh->has_tag(VERT, "metric")) {
      return measure_edges_tmpl<MetricEdgeLengths<2>>(mesh, a2e);
    }
  }
  osh_fail("measure_edges(): no size field exists!\n");
  NORETURN(Reals());
}

Reals measure_edges_real(Mesh* mesh) {
  return measure_edges_real(mesh, LOs(mesh->nedges(), 0, 1));
}

Reals measure_edges_metric(Mesh* mesh) {
  return measure_edges_metric(mesh, LOs(mesh->nedges(), 0, 1));
}

Reals find_identity_size(Mesh* mesh) {
  CHECK(mesh->owners_have_all_upward(VERT));
  auto lens = measure_edges_real(mesh);
  auto v2e = mesh->ask_up(VERT, EDGE);
  auto nve = v2e.a2ab.last();
  auto weights = Reals(nve, 1.0);
  auto own_isos = graph_weighted_average(v2e, weights, lens, 1);
  auto synced_isos = mesh->sync_array(VERT, own_isos, 1);
  return synced_isos;
}

template <Int dim>
static Reals measure_elements_real_tmpl(Mesh* mesh) {
  RealElementSizes measurer(mesh);
  auto ev2v = mesh->ask_verts_of(dim);
  auto ne = mesh->nelems();
  Write<Real> sizes(ne);
  auto f = LAMBDA(LO e) {
    auto v = gather_verts<dim + 1>(ev2v, e);
    sizes[e] = measurer.measure(v);
  };
  parallel_for(ne, f);
  return sizes;
}

Reals measure_elements_real(Mesh* mesh) {
  if (mesh->dim() == 3) {
    return measure_elements_real_tmpl<3>(mesh);
  } else {
    CHECK(mesh->dim() == 2);
    return measure_elements_real_tmpl<2>(mesh);
  }
}

template <Int dim>
static Reals find_identity_metric_tmpl(Mesh* mesh) {
  CHECK(dim == mesh->dim());
  CHECK(mesh->owners_have_all_upward(VERT));
  auto coords = mesh->coords();
  auto ev2v = mesh->ask_verts_of(dim);
  auto elem_metrics_w = Write<Real>(mesh->nelems() * symm_dofs(dim));
  auto f0 = LAMBDA(LO e) {
    auto v = gather_verts<dim + 1>(ev2v, e);
    auto p = gather_vectors<dim + 1, dim>(coords, v);
    auto m = element_identity_metric(p);
    set_symm(elem_metrics_w, e, m);
  };
  parallel_for(mesh->nelems(), f0);
  auto elem_metrics = Reals(elem_metrics_w);
  auto elem_sizes = measure_elements_real(mesh);
  auto v2e = mesh->ask_up(VERT, dim);
  auto v2ve = v2e.a2ab;
  auto ve2e = v2e.ab2b;
  auto vert_metrics_w = Write<Real>(mesh->nverts() * symm_dofs(dim));
  auto f1 = LAMBDA(LO v) {
    Real ess = 0.0;
    auto iems = zero_matrix<dim, dim>();
    for (auto ve = v2ve[v]; ve < v2ve[v + 1]; ++ve) {
      auto e = ve2e[ve];
      auto em = get_symm<dim>(elem_metrics, e);
      auto es = elem_sizes[e];
      auto iem = invert(em);
      iems = iems + (iem * es);
      ess += es;
    }
    auto vm = invert(iems / ess);
    set_symm(vert_metrics_w, v, vm);
  };
  parallel_for(mesh->nverts(), f1);
  return vert_metrics_w;
}

Reals find_identity_metric(Mesh* mesh) {
  if (mesh->dim() == 3) {
    return find_identity_metric_tmpl<3>(mesh);
  } else {
    CHECK(mesh->dim() == 2);
    return find_identity_metric_tmpl<2>(mesh);
  }
}

/* The algorithms below are for scaling a size field such that
 * adapting based on that size field will result in a certain specified
 * number of elements.
 * The formulas are based on Section 2.7 of:
 * Pain, C. C., et al.
 * "Tetrahedral mesh optimisation and adaptivity for
 *  steady-state and transient finite element calculations."
 * Computer Methods in Applied Mechanics and Engineering
 * 190.29 (2001): 3771-3796.
 *
 * We suspect that Pain et al.'s $\theta$ correction factor is needed
 * because the tetrahedra in a real mesh are not equilateral, and they
 * do use equilateral volume in their formula.
 * To correct for this, we will use information available in our
 * mean ratio shape measure.
 * In particular, the following should hold:
 *
 * $\eta^\frac{3}{2} = \frac{1}{V_{eq}} \frac{V_e}{l_{RMS}^3}$
 *
 * Where $\eta$ is the mean ratio (see quality.hpp),
 * $V_{eq}$ is the volume of an equilateral tetrahedron with unit
 * edge lengths, $V_e$ is the volume of the element being measured,
 * and $l_{RMS}$ is its root-mean-squared edge length.
 * Loosely speaking, the mean ratio to some power is the volume
 * that this tetrahedron would have if we scaled it until its
 * root mean squared edge length was one.
 * This happens to be exactly the correction factor we need in this case.
 *
 * When we need to scale desired edge lengths, we
 * use the inverse of the square root of the $\beta$ factor.
 * (a metric tensor eigenvalue is $(1/h^2)$, where $h$ is desired length).
*/

struct IsoDilation {
  Reals e2h;
  IsoDilation(Mesh* mesh, Reals v2h) {
    auto e2e = LOs(mesh->nelems(), 0, 1);
    e2h = average_field(mesh, mesh->dim(), e2e, 1, v2h);
  }
  template <Int dim>
  DEVICE Real get_dilation(LO e) const {
    return 1. / raise(e2h[e], dim);
  }
};

struct AnisoDilation {
  Reals e2m;
  AnisoDilation(Mesh* mesh, Reals v2m) {
    auto e2e = LOs(mesh->nelems(), 0, 1);
    e2m = average_metric(mesh, mesh->dim(), e2e, v2m);
  }
  template <Int dim>
  DEVICE Real get_dilation(LO e) const {
    return sqrt(determinant(get_symm<dim>(e2m, e)));
  }
};

template <Int dim, typename Dilation>
static Reals perfect_metric_volumes_tmpl(Mesh* mesh, Reals v2sf) {
  auto dilation_measurer = Dilation(mesh, v2sf);
  auto e2q = mesh->ask_qualities();
  auto ev2v = mesh->ask_verts_of(dim);
  auto real_vol_measurer = RealElementSizes(mesh);
  auto out_w = Write<Real>(mesh->nelems());
  auto f = LAMBDA(LO e) {
    auto evv2v = gather_verts<dim + 1>(ev2v, e);
    auto dilation = dilation_measurer.template get_dilation<dim>(e);
    auto volume = real_vol_measurer.measure(evv2v);
    auto quality = e2q[e];
    auto quality_correction = 1. / sqrt(raise(quality, dim));
    out_w[e] = volume * dilation * quality_correction;
  };
  parallel_for(mesh->nelems(), f);
  return Reals(out_w);
}

Reals perfect_size_volumes(Mesh* mesh, Reals v2h) {
  if (mesh->dim() == 3) return perfect_metric_volumes_tmpl<3,IsoDilation>(mesh, v2h);
  if (mesh->dim() == 2) return perfect_metric_volumes_tmpl<2,IsoDilation>(mesh, v2h);
  NORETURN(Reals());
}

Reals perfect_metric_volumes(Mesh* mesh, Reals v2m) {
  if (mesh->dim() == 3) return perfect_metric_volumes_tmpl<3,AnisoDilation>(mesh, v2m);
  if (mesh->dim() == 2) return perfect_metric_volumes_tmpl<2,AnisoDilation>(mesh, v2m);
  NORETURN(Reals());
}

Real volume_scalar_for_nelems(Mesh* mesh, Real volume_sum, Real target_nelems) {
  Real unit_edge_volume = 0;
  if (mesh->dim() == 3) unit_edge_volume = 1. / sqrt(72.);
  if (mesh->dim() == 2) unit_edge_volume = sqrt(3.) / 4.;
  return target_nelems * unit_edge_volume / volume_sum;
}

Reals scale_size_for_nelems(Mesh* mesh, Reals v2h, Real target_nelems) {
  auto vols = perfect_size_volumes(mesh, v2h);
  auto vol_sum = repro_sum_owned(mesh, mesh->dim(), vols);
  auto vol_scal = volume_scalar_for_nelems(mesh, vol_sum, target_nelems);
  Real h_scal = 0;
  if (mesh->dim() == 3) h_scal = 1. / cbrt(vol_scal);
  if (mesh->dim() == 2) h_scal = 1. / sqrt(vol_scal);
  return multiply_each_by(h_scal, v2h);
}

Reals scale_metric_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems) {
  auto vols = perfect_metric_volumes(mesh, v2m);
  auto vol_sum = repro_sum_owned(mesh, mesh->dim(), vols);
  auto vol_scal = volume_scalar_for_nelems(mesh, vol_sum, target_nelems);
  Real m_scal = 0;
  if (mesh->dim() == 3) m_scal = cbrt(square(vol_scal));
  if (mesh->dim() == 2) m_scal = vol_scal;
  return multiply_each_by(m_scal, v2m);
}

}  // end namespace osh
