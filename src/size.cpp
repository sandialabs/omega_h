#include "size.hpp"

#include <iostream>

#include "array.hpp"
#include "eigen.hpp"
#include "graph.hpp"
#include "loop.hpp"
#include "project.hpp"
#include "quality.hpp"

namespace Omega_h {

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
  Omega_h_fail("measure_edges(): no size field exists!\n");
  NORETURN(Reals());
}

Reals measure_edges_real(Mesh* mesh) {
  return measure_edges_real(mesh, LOs(mesh->nedges(), 0, 1));
}

Reals measure_edges_metric(Mesh* mesh) {
  return measure_edges_metric(mesh, LOs(mesh->nedges(), 0, 1));
}

template <Int dim>
static Reals measure_elements_real_tmpl(Mesh* mesh) {
  RealElementSizes measurer(mesh);
  auto ev2v = mesh->ask_elem_verts();
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
  if (mesh->dim() == 3) return measure_elements_real_tmpl<3>(mesh);
  if (mesh->dim() == 2) return measure_elements_real_tmpl<2>(mesh);
  NORETURN(Reals());
}

template <Int dim>
static Reals element_implied_sizes_dim(Mesh* mesh) {
  auto coords = mesh->coords();
  auto ev2v = mesh->ask_elem_verts();
  auto out = Write<Real>(mesh->nelems());
  auto f = LAMBDA(LO e) {
    auto v = gather_verts<dim + 1>(ev2v, e);
    auto p = gather_vectors<dim + 1, dim>(coords, v);
    auto h = element_implied_size(p);
    out[e] = h;
  };
  parallel_for(mesh->nelems(), f);
  return out;
}

static Reals element_implied_sizes(Mesh* mesh) {
  if (mesh->dim() == 3) return element_implied_sizes_dim<3>(mesh);
  if (mesh->dim() == 2) return element_implied_sizes_dim<2>(mesh);
  NORETURN(Reals());
}

Reals find_implied_size(Mesh* mesh) {
  auto e_h = element_implied_sizes(mesh);
  auto e_linear = linearize_isos(e_h);
  auto v_linear = project_by_average(mesh, e_linear);
  auto v_h = delinearize_isos(v_linear);
  return v_h;
}

template <Int dim>
static Reals element_implied_metrics_dim(Mesh* mesh) {
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nelems() * symm_dofs(dim));
  auto f = LAMBDA(LO e) {
    auto v = gather_verts<dim + 1>(ev2v, e);
    auto p = gather_vectors<dim + 1, dim>(coords, v);
    auto m = element_implied_metric(p);
    set_symm(out, e, m);
  };
  parallel_for(mesh->nelems(), f);
  return out;
}

static Reals element_implied_metrics(Mesh* mesh) {
  if (mesh->dim() == 3) return element_implied_metrics_dim<3>(mesh);
  if (mesh->dim() == 2) return element_implied_metrics_dim<2>(mesh);
  NORETURN(Reals());
}

Reals find_implied_metric(Mesh* mesh) {
  auto e_metric = element_implied_metrics(mesh);
  auto e_linear = linearize_metrics(mesh->dim(), e_metric);
  auto v_linear = project_by_average(mesh, e_linear);
  auto v_metric = delinearize_metrics(mesh->dim(), v_linear);
  return v_metric;
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
 * In fact, if we use such a correction factor, most terms cancel
 * out and we are left with the following weights for existing elements:
 *
 * isotropic case: $\frac{l_{RMS}^3}{h^3}$
 * anisotropic case: $\frac{l_{M,RMS}^3}$
 *
 * Where $l_{M,RMS}$ is the root mean squared value of $v_i^T \mathcal{M} v_i$,
 * for edge vector $v_i$, i.e. the root mean squared metric edge length.
 * In both cases, the (more accurate) scaling factor is simply
 * the root mean squared edge length in metric space
 * to the power of the element dimension.
 *
 * When we need to scale desired edge lengths, we
 * use the inverse of the square root of the $\beta$ factor.
 * (a metric tensor eigenvalue is $(1/h^2)$, where $h$ is desired length).
*/

template <typename EdgeVectors>
INLINE Real mean_squared_iso_length(EdgeVectors edge_vectors, Real h) {
  auto nedges = EdgeVectors::size;
  Real msl = 0;
  for (Int i = 0; i < nedges; ++i) {
    msl += norm_squared(edge_vectors[i]) / square(h);
  }
  return msl / nedges;
}

struct MeanSquaredIsoLength {
  Reals e2h;
  MeanSquaredIsoLength(Mesh* mesh, Reals v2h) {
    auto e2e = LOs(mesh->nelems(), 0, 1);
    e2h = average_field(mesh, mesh->dim(), e2e, 1, v2h);
  }
  template <Int dim, typename EdgeVectors>
  DEVICE Real get(LO e, EdgeVectors edge_vectors) const {
    auto h = e2h[e];
    return mean_squared_iso_length(edge_vectors, h);
  }
};

struct MeanSquaredMetricLength {
  Reals e2m;
  MeanSquaredMetricLength(Mesh* mesh, Reals v2m) {
    auto e2e = LOs(mesh->nelems(), 0, 1);
    /* TODO: if this is the only use of get_maxdet_metrics, we can
     * omit the dummy e2e map altogether
     */
    e2m = get_maxdet_metrics(mesh, e2e, v2m);
  }
  template <Int dim, typename EdgeVectors>
  DEVICE Real get(LO e, EdgeVectors edge_vectors) const {
    auto m = get_symm<dim>(e2m, e);
    return mean_squared_metric_length(edge_vectors, m);
  }
};

template <Int dim, typename MeanSquaredLength>
static Reals expected_elems_per_elem_tmpl(Mesh* mesh, Reals v2sf) {
  auto msl_obj = MeanSquaredLength(mesh, v2sf);
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto out_w = Write<Real>(mesh->nelems());
  auto f = LAMBDA(LO e) {
    auto eev2v = gather_verts<dim + 1>(ev2v, e);
    auto eev2x = gather_vectors<dim + 1, dim>(coords, eev2v);
    auto basis = simplex_basis<dim, dim>(eev2x);
    auto edge_vectors = element_edge_vectors(eev2x, basis);
    auto msl = msl_obj.template get<dim>(e, edge_vectors);
    out_w[e] = power<dim, 2>(msl);
  };
  parallel_for(mesh->nelems(), f);
  return Reals(out_w);
}

Reals expected_elems_per_elem_iso(Mesh* mesh, Reals v2h) {
  if (mesh->dim() == 3) {
    return expected_elems_per_elem_tmpl<3, MeanSquaredIsoLength>(mesh, v2h);
  }
  if (mesh->dim() == 2) {
    return expected_elems_per_elem_tmpl<2, MeanSquaredIsoLength>(mesh, v2h);
  }
  NORETURN(Reals());
}

Reals expected_elems_per_elem_metric(Mesh* mesh, Reals v2m) {
  if (mesh->dim() == 3) {
    return expected_elems_per_elem_tmpl<3, MeanSquaredMetricLength>(mesh, v2m);
  }
  if (mesh->dim() == 2) {
    return expected_elems_per_elem_tmpl<2, MeanSquaredMetricLength>(mesh, v2m);
  }
  NORETURN(Reals());
}

Real size_scalar_for_nelems(Mesh* mesh, Reals v2h, Real target_nelems) {
  auto elems_per_elem = expected_elems_per_elem_iso(mesh, v2h);
  auto elems = repro_sum_owned(mesh, mesh->dim(), elems_per_elem);
  auto size_scal = target_nelems / elems;
  Real h_scal = 0;
  if (mesh->dim() == 3) h_scal = 1. / cbrt(size_scal);
  if (mesh->dim() == 2) h_scal = 1. / sqrt(size_scal);
  return h_scal;
}

Real metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems) {
  auto elems_per_elem = expected_elems_per_elem_metric(mesh, v2m);
  auto elems = repro_sum_owned(mesh, mesh->dim(), elems_per_elem);
  auto size_scal = target_nelems / elems;
  Real m_scal = 0;
  if (mesh->dim() == 3) m_scal = cbrt(square(size_scal));
  if (mesh->dim() == 2) m_scal = size_scal;
  return m_scal;
}

Reals get_midedge_isos(Mesh* mesh, LOs a2e, Reals v2h) {
  auto na = a2e.size();
  Write<Real> out(na);
  auto ev2v = mesh->ask_verts_of(EDGE);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<2>(ev2v, e);
    auto hs = gather_scalars<2>(v2h, v);
    auto h = interpolate_metric(hs[0], hs[1], 1.0 / 2.0);
    out[a] = h;
  };
  parallel_for(na, f);
  return out;
}

Reals interpolate_between_isos(Reals a, Reals b, Real t) {
  auto log_a = linearize_isos(a);
  auto log_b = linearize_isos(b);
  auto log_c = interpolate_between(a, b, t);
  return delinearize_isos(log_c);
}

Reals linearize_isos(Reals isos) {
  auto n = isos.size();
  auto out = Write<Real>(n);
  auto f = LAMBDA(LO i) { out[i] = linearize_metric(isos[i]); };
  parallel_for(n, f);
  return out;
}

Reals delinearize_isos(Reals log_isos) {
  auto n = log_isos.size();
  auto out = Write<Real>(n);
  auto f = LAMBDA(LO i) { out[i] = delinearize_metric(log_isos[i]); };
  parallel_for(n, f);
  return out;
}

}  // end namespace Omega_h
