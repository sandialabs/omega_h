#include "Omega_h_metric.hpp"

#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_recover.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_simplex.hpp"
#include "host_few.hpp"
#include "surface.hpp"

namespace Omega_h {

Int get_metric_dim(Int ncomps) {
  for (Int i = 1; i <= 3; ++i)
    if (ncomps == symm_ncomps(i)) return i;
  NORETURN(Int());
}

Int get_metrics_dim(LO nmetrics, Reals metrics) {
  auto ncomps = divide_no_remainder(metrics.size(), nmetrics);
  return get_metric_dim(ncomps);
}

Int get_metric_dim(Mesh* mesh) {
  auto ncomps = mesh->get_tagbase(VERT, "metric")->ncomps();
  return get_metric_dim(ncomps);
}

template <Int dim>
static Reals clamp_metrics_dim(
    LO nmetrics, Reals metrics, Real h_min, Real h_max) {
  auto out = Write<Real>(nmetrics * symm_ncomps(dim));
  auto f = LAMBDA(LO i) {
    auto m = get_symm<dim>(metrics, i);
    m = clamp_metric(m, h_min, h_max);
    set_symm(out, i, m);
  };
  parallel_for(nmetrics, f);
  return out;
}

Reals clamp_metrics(LO nmetrics, Reals metrics, Real h_min, Real h_max) {
  auto dim = get_metrics_dim(nmetrics, metrics);
  if (dim == 3) return clamp_metrics_dim<3>(nmetrics, metrics, h_min, h_max);
  if (dim == 2) return clamp_metrics_dim<2>(nmetrics, metrics, h_min, h_max);
  if (dim == 1) return clamp_metrics_dim<1>(nmetrics, metrics, h_min, h_max);
  NORETURN(Reals());
}

template <Int mdim, Int edim>
static Reals mident_metrics_tmpl(Mesh* mesh, LOs a2e, Reals v2m) {
  auto na = a2e.size();
  Write<Real> out(na * symm_ncomps(mdim));
  auto ev2v = mesh->ask_verts_of(edim);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<edim + 1>(ev2v, e);
    auto ms = gather_symms<edim + 1, mdim>(v2m, v);
    auto m = average_metric(ms);
    set_symm(out, a, m);
  };
  parallel_for(na, f);
  return out;
}

Reals get_mident_metrics(Mesh* mesh, Int ent_dim, LOs entities, Reals v2m) {
  auto metrics_dim = get_metrics_dim(mesh->nverts(), v2m);
  if (metrics_dim == 3 && ent_dim == 3) {
    return mident_metrics_tmpl<3, 3>(mesh, entities, v2m);
  }
  if (metrics_dim == 3 && ent_dim == 1) {
    return mident_metrics_tmpl<3, 1>(mesh, entities, v2m);
  }
  if (metrics_dim == 2 && ent_dim == 2) {
    return mident_metrics_tmpl<2, 2>(mesh, entities, v2m);
  }
  if (metrics_dim == 2 && ent_dim == 1) {
    return mident_metrics_tmpl<2, 1>(mesh, entities, v2m);
  }
  if (metrics_dim == 1 && ent_dim == 3) {
    return mident_metrics_tmpl<1, 3>(mesh, entities, v2m);
  }
  if (metrics_dim == 1 && ent_dim == 2) {
    return mident_metrics_tmpl<1, 2>(mesh, entities, v2m);
  }
  if (metrics_dim == 1 && ent_dim == 1) {
    return mident_metrics_tmpl<1, 1>(mesh, entities, v2m);
  }
  NORETURN(Reals());
}

Reals interpolate_between_metrics(LO nmetrics, Reals a, Reals b, Real t) {
  auto log_a = linearize_metrics(nmetrics, a);
  auto log_b = linearize_metrics(nmetrics, b);
  auto log_c = interpolate_between(log_a, log_b, t);
  return delinearize_metrics(nmetrics, log_c);
}

template <Int dim>
Reals linearize_metrics_dim(Reals metrics) {
  auto n = metrics.size() / symm_ncomps(dim);
  auto out = Write<Real>(n * symm_ncomps(dim));
  auto f = LAMBDA(LO i) {
    set_symm(out, i, linearize_metric(get_symm<dim>(metrics, i)));
  };
  parallel_for(n, f);
  return out;
}

template <Int dim>
Reals delinearize_metrics_dim(Reals lms) {
  auto n = lms.size() / symm_ncomps(dim);
  auto out = Write<Real>(n * symm_ncomps(dim));
  auto f = LAMBDA(LO i) {
    set_symm(out, i, delinearize_metric(get_symm<dim>(lms, i)));
  };
  parallel_for(n, f);
  return out;
}

Reals linearize_metrics(LO nmetrics, Reals metrics) {
  auto dim = get_metrics_dim(nmetrics, metrics);
  if (dim == 3) return linearize_metrics_dim<3>(metrics);
  if (dim == 2) return linearize_metrics_dim<2>(metrics);
  if (dim == 1) return linearize_metrics_dim<1>(metrics);
  NORETURN(Reals());
}

Reals delinearize_metrics(LO nmetrics, Reals linear_metrics) {
  auto dim = get_metrics_dim(nmetrics, linear_metrics);
  if (dim == 3) return delinearize_metrics_dim<3>(linear_metrics);
  if (dim == 2) return delinearize_metrics_dim<2>(linear_metrics);
  if (dim == 1) return delinearize_metrics_dim<1>(linear_metrics);
  NORETURN(Reals());
}

template <Int dim>
static HostFew<Reals, dim> axes_from_metrics_dim(Reals metrics) {
  auto n = divide_no_remainder(metrics.size(), symm_ncomps(dim));
  HostFew<Write<Real>, dim> w;
  for (Int i = 0; i < dim; ++i) w[i] = Write<Real>(n * dim);
  auto f = LAMBDA(LO i) {
    auto md = decompose_metric(get_symm<dim>(metrics, i));
    for (Int j = 0; j < dim; ++j) set_vector(w[j], i, md.q[j] * md.l[j]);
  };
  parallel_for(n, f);
  HostFew<Reals, dim> r;
  for (Int i = 0; i < dim; ++i) r[i] = Reals(w[i]);
  return r;
}

template <Int dim>
static void axes_from_metric_field_dim(Mesh* mesh,
    std::string const& metric_name, std::string const& output_prefix) {
  auto metrics = mesh->get_array<Real>(VERT, metric_name);
  auto axes = axes_from_metrics_dim<dim>(metrics);
  for (Int i = 0; i < dim; ++i) {
    mesh->add_tag(VERT, output_prefix + '_' + std::to_string(i), dim, axes[i]);
  }
}

void axes_from_metric_field(Mesh* mesh, std::string const& metric_name,
    std::string const& axis_prefix) {
  if (mesh->dim() == 3) {
    axes_from_metric_field_dim<3>(mesh, metric_name, axis_prefix);
    return;
  }
  if (mesh->dim() == 2) {
    axes_from_metric_field_dim<2>(mesh, metric_name, axis_prefix);
    return;
  }
  NORETURN();
}

/* gradation limiting code: */

template <Int mesh_dim, Int metric_dim>
static Reals limit_gradation_once_tmpl(
    Mesh* mesh, Reals values, Real max_rate) {
  auto v2v = mesh->ask_star(VERT);
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(metric_dim));
  auto f = LAMBDA(LO v) {
    auto m = get_symm<metric_dim>(values, v);
    auto x = get_vector<mesh_dim>(coords, v);
    for (auto vv = v2v.a2ab[v]; vv < v2v.a2ab[v + 1]; ++vv) {
      auto av = v2v.ab2b[vv];
      auto am = get_symm<metric_dim>(values, av);
      auto ax = get_vector<mesh_dim>(coords, av);
      auto vec = ax - x;
      auto metric_dist = metric_length(m, vec);
      auto decomp = decompose_metric(m);
      decomp.l = decomp.l * (1.0 + metric_dist * max_rate);
      auto limiter = compose_metric(decomp.q, decomp.l);
      auto limited = intersect_metrics(m, limiter);
      m = limited;
    }
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  values = Reals(out);
  values = mesh->sync_array(VERT, values, symm_ncomps(metric_dim));
  return values;
}

static Reals limit_gradation_once(Mesh* mesh, Reals values, Real max_rate) {
  auto metric_dim = get_metrics_dim(mesh->nverts(), values);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return limit_gradation_once_tmpl<3, 3>(mesh, values, max_rate);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return limit_gradation_once_tmpl<2, 2>(mesh, values, max_rate);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return limit_gradation_once_tmpl<3, 1>(mesh, values, max_rate);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return limit_gradation_once_tmpl<2, 1>(mesh, values, max_rate);
  }
  NORETURN(Reals());
}

Reals limit_metric_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol) {
  CHECK(mesh->owners_have_all_upward(VERT));
  CHECK(max_rate > 0.0);
  auto comm = mesh->comm();
  Reals values2 = values;
  Int i = 0;
  do {
    values = values2;
    values2 = limit_gradation_once(mesh, values, max_rate);
    ++i;
    if (can_print(mesh) && i > 40) {
      std::cout << "warning: gradation limiting is up to step " << i << '\n';
    }
  } while (!comm->reduce_and(are_close(values, values2, tol)));
  if (can_print(mesh)) {
    std::cout << "limited gradation in " << i << " steps\n";
  }
  return values2;
}

Reals project_metrics(Mesh* mesh, Reals e2m) {
  auto e_linear = linearize_metrics(mesh->nelems(), e2m);
  auto v_linear = project_by_average(mesh, e_linear);
  return delinearize_metrics(mesh->nverts(), v_linear);
}

Reals smooth_metric_once(Mesh* mesh, Reals v2m) {
  auto e2e = LOs(mesh->nelems(), 0, 1);
  return project_metrics(mesh, get_mident_metrics(mesh, mesh->dim(), e2e, v2m));
}

template <Int dim>
static Reals element_implied_isos_dim(Mesh* mesh) {
  auto coords = mesh->coords();
  auto ev2v = mesh->ask_elem_verts();
  auto out = Write<Real>(mesh->nelems());
  auto f = LAMBDA(LO e) {
    auto v = gather_verts<dim + 1>(ev2v, e);
    auto p = gather_vectors<dim + 1, dim>(coords, v);
    auto h = element_implied_length(p);
    out[e] = metric_eigenvalue_from_length(h);
  };
  parallel_for(mesh->nelems(), f);
  return out;
}

static Reals element_implied_isos(Mesh* mesh) {
  if (mesh->dim() == 3) return element_implied_isos_dim<3>(mesh);
  if (mesh->dim() == 2) return element_implied_isos_dim<2>(mesh);
  NORETURN(Reals());
}

Reals find_implied_isos(Mesh* mesh) {
  return project_metrics(mesh, element_implied_isos(mesh));
}

template <Int dim>
static Reals element_implied_metrics_dim(Mesh* mesh) {
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nelems() * symm_ncomps(dim));
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
  return project_metrics(mesh, element_implied_metrics(mesh));
}

/* A Hessian-based anisotropic size field, from
 * Alauzet's tech report:
 *
 * F. Alauzet, P.J. Frey, Estimateur d'erreur geometrique
 * et metriques anisotropes pour l'adaptation de maillage.
 * Partie I: aspects theoriques,
 * RR-4759, INRIA Rocquencourt, 2003.
 */

template <Int dim>
static INLINE Matrix<dim, dim> metric_from_hessian(
    Matrix<dim, dim> hessian, Real eps) {
  auto ed = decompose_eigen(hessian);
  auto r = ed.q;
  auto l = ed.l;
  constexpr auto c_num = square(dim);
  constexpr auto c_denom = 2 * square(dim + 1);
  decltype(l) tilde_l;
  for (Int i = 0; i < dim; ++i) {
    tilde_l[i] = (c_num * fabs(l[i])) / (c_denom * eps);
  }
  return compose_eigen(r, tilde_l);
}

template <Int dim>
static Reals metric_from_hessians_dim(Reals hessians, Real eps) {
  auto ncomps = symm_ncomps(dim);
  auto n = divide_no_remainder(hessians.size(), ncomps);
  auto out = Write<Real>(n * ncomps);
  auto f = LAMBDA(LO i) {
    auto hess = get_symm<dim>(hessians, i);
    auto m = metric_from_hessian(hess, eps);
    set_symm(out, i, m);
  };
  parallel_for(n, f);
  return out;
}

Reals metric_from_hessians(Int dim, Reals hessians, Real eps) {
  CHECK(eps > 0.0);
  if (dim == 3) return metric_from_hessians_dim<3>(hessians, eps);
  if (dim == 2) return metric_from_hessians_dim<2>(hessians, eps);
  NORETURN(Reals());
}

Reals get_curvature_isos(Mesh* mesh, Real segment_angle) {
  auto vert_curvatures = get_vert_curvatures(mesh);
  auto out = Write<Real>(mesh->nverts());
  auto f = LAMBDA(LO v) {
    auto curvature = vert_curvatures[v];
    auto l = square(curvature / segment_angle);
    out[v] = l;
  };
  parallel_for(mesh->nverts(), f);
  return out;
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

struct MeanSquaredMetricLength {
  Reals e2m;
  MeanSquaredMetricLength(Mesh* mesh, Reals v2m) {
    auto e2e = LOs(mesh->nelems(), 0, 1);
    e2m = get_mident_metrics(mesh, mesh->dim(), e2e, v2m);
  }
  template <Int metric_dim, typename EdgeVectors>
  DEVICE Real get(LO e, EdgeVectors edge_vectors) const {
    auto m = get_symm<metric_dim>(e2m, e);
    return mean_squared_metric_length(edge_vectors, m);
  }
};

template <Int mesh_dim, Int metric_dim>
static Reals expected_elems_per_elem_tmpl(Mesh* mesh, Reals v2m) {
  auto msl_obj = MeanSquaredMetricLength(mesh, v2m);
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto out_w = Write<Real>(mesh->nelems());
  auto f = LAMBDA(LO e) {
    auto eev2v = gather_verts<mesh_dim + 1>(ev2v, e);
    auto eev2x = gather_vectors<mesh_dim + 1, mesh_dim>(coords, eev2v);
    auto basis = simplex_basis<mesh_dim, mesh_dim>(eev2x);
    auto edge_vectors = element_edge_vectors(eev2x, basis);
    auto msl = msl_obj.template get<metric_dim>(e, edge_vectors);
    out_w[e] = power<mesh_dim, 2>(msl);
  };
  parallel_for(mesh->nelems(), f);
  return Reals(out_w);
}

Reals expected_elems_per_elem(Mesh* mesh, Reals v2m) {
  auto metric_dim = get_metrics_dim(mesh->nverts(), v2m);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return expected_elems_per_elem_tmpl<3, 3>(mesh, v2m);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return expected_elems_per_elem_tmpl<2, 2>(mesh, v2m);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return expected_elems_per_elem_tmpl<3, 1>(mesh, v2m);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return expected_elems_per_elem_tmpl<2, 1>(mesh, v2m);
  }
  NORETURN(Reals());
}

Real metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems) {
  auto elems_per_elem = expected_elems_per_elem(mesh, v2m);
  auto elems = repro_sum_owned(mesh, mesh->dim(), elems_per_elem);
  auto size_scal = target_nelems / elems;
  auto metric_dim = get_metrics_dim(mesh->nverts(), v2m);
  auto metric_scal = power(size_scal, 2, metric_dim);
  return metric_scal;
}

}  // end namespace Omega_h
