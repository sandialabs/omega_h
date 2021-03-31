#include "Omega_h_metric.hpp"

#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_confined.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_host_few.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_metric_intersect.hpp"
#include "Omega_h_recover.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_simplex.hpp"
#include "Omega_h_surface.hpp"

namespace Omega_h {

Int get_metric_dim(Int ncomps) {
  for (Int i = 1; i <= 3; ++i)
    if (ncomps == symm_ncomps(i)) return i;
  OMEGA_H_NORETURN(Int());
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
Reals clamp_metrics_dim(
    LO nmetrics, Reals metrics, Real h_min, Real h_max) {
  auto out = Write<Real>(nmetrics * symm_ncomps(dim));
  auto functor = OMEGA_H_LAMBDA(LO i) {
    auto m = get_symm<dim>(metrics, i);
    m = clamp_metric(m, h_min, h_max);
    set_symm(out, i, m);
  };
  parallel_for("clamp_metrics", nmetrics, std::move(functor));
  return out;
}

Reals clamp_metrics(LO nmetrics, Reals metrics, Real h_min, Real h_max) {
  auto dim = get_metrics_dim(nmetrics, metrics);
  if (dim == 3) return clamp_metrics_dim<3>(nmetrics, metrics, h_min, h_max);
  if (dim == 2) return clamp_metrics_dim<2>(nmetrics, metrics, h_min, h_max);
  if (dim == 1) return clamp_metrics_dim<1>(nmetrics, metrics, h_min, h_max);
  OMEGA_H_NORETURN(Reals());
}

template <Int mdim, Int edim>
 Reals get_mident_metrics_tmpl(
    Mesh* mesh, LOs a2e, Reals v2m, bool has_degen) {
  auto na = a2e.size();
  Write<Real> out(na * symm_ncomps(mdim));
  auto ev2v = mesh->ask_verts_of(edim);
  auto f = OMEGA_H_LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<edim + 1>(ev2v, e);
    auto ms = gather_symms<edim + 1, mdim>(v2m, v);
    auto m = average_metric(ms, has_degen);
    set_symm(out, a, m);
  };
  parallel_for(na, f, "get_mident_metrics");
  return out;
}

Reals get_mident_metrics(
    Mesh* mesh, Int ent_dim, LOs entities, Reals v2m, bool has_degen) {
  if (entities.size() == 0) return Reals({});
  auto metrics_dim = get_metrics_dim(mesh->nverts(), v2m);
  if (metrics_dim == 3 && ent_dim == 3) {
    return get_mident_metrics_tmpl<3, 3>(mesh, entities, v2m, has_degen);
  }
  if (metrics_dim == 3 && ent_dim == 1) {
    return get_mident_metrics_tmpl<3, 1>(mesh, entities, v2m, has_degen);
  }
  if (metrics_dim == 2 && ent_dim == 2) {
    return get_mident_metrics_tmpl<2, 2>(mesh, entities, v2m, has_degen);
  }
  if (metrics_dim == 2 && ent_dim == 1) {
    return get_mident_metrics_tmpl<2, 1>(mesh, entities, v2m, has_degen);
  }
  if (metrics_dim == 1 && ent_dim == 3) {
    return get_mident_metrics_tmpl<1, 3>(mesh, entities, v2m, has_degen);
  }
  if (metrics_dim == 1 && ent_dim == 2) {
    return get_mident_metrics_tmpl<1, 2>(mesh, entities, v2m, has_degen);
  }
  if (metrics_dim == 1 && ent_dim == 1) {
    return get_mident_metrics_tmpl<1, 1>(mesh, entities, v2m, has_degen);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals get_mident_metrics(Mesh* mesh, Int ent_dim, Reals v2m, bool has_degen) {
  LOs e2e(mesh->nents(ent_dim), 0, 1);
  return get_mident_metrics(mesh, ent_dim, e2e, v2m, has_degen);
}

Reals interpolate_between_metrics(LO nmetrics, Reals a, Reals b, Real t) {
  auto log_a = linearize_metrics(nmetrics, a);
  auto log_b = linearize_metrics(nmetrics, b);
  auto log_c = interpolate_between(log_a, log_b, t);
  return delinearize_metrics(nmetrics, log_c);
}

template <Int dim>
Reals linearize_metrics_dim(Reals metrics) {
  auto n = divide_no_remainder(metrics.size(), symm_ncomps(dim));
  auto out = Write<Real>(n * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO i) {
    set_symm(out, i, linearize_metric(get_symm<dim>(metrics, i)));
  };
  parallel_for(n, f, "linearize_metrics");
  return out;
}

template <Int dim>
Reals delinearize_metrics_dim(Reals lms) {
  auto n = divide_no_remainder(lms.size(), symm_ncomps(dim));
  auto out = Write<Real>(n * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO i) {
    set_symm(out, i, delinearize_metric(get_symm<dim>(lms, i)));
  };
  parallel_for(n, f, "delinearize_metrics");
  return out;
}

Reals linearize_metrics(LO nmetrics, Reals metrics) {
  if (nmetrics == 0) return metrics;
  auto dim = get_metrics_dim(nmetrics, metrics);
  if (dim == 3) return linearize_metrics_dim<3>(metrics);
  if (dim == 2) return linearize_metrics_dim<2>(metrics);
  if (dim == 1) return linearize_metrics_dim<1>(metrics);
  OMEGA_H_NORETURN(Reals());
}

Reals delinearize_metrics(LO nmetrics, Reals linear_metrics) {
  if (nmetrics == 0) return linear_metrics;
  auto dim = get_metrics_dim(nmetrics, linear_metrics);
  if (dim == 3) return delinearize_metrics_dim<3>(linear_metrics);
  if (dim == 2) return delinearize_metrics_dim<2>(linear_metrics);
  if (dim == 1) return delinearize_metrics_dim<1>(linear_metrics);
  OMEGA_H_NORETURN(Reals());
}

/* gradation limiting code: */

template <Int mesh_dim, Int metric_dim>
Reals limit_gradation_once_tmpl(
    Mesh* mesh, Reals values, Real max_rate) {
  auto v2v = mesh->ask_star(VERT);
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(metric_dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto m = get_symm<metric_dim>(values, v);
    auto x = get_vector<mesh_dim>(coords, v);
    for (auto vv = v2v.a2ab[v]; vv < v2v.a2ab[v + 1]; ++vv) {
      auto av = v2v.ab2b[vv];
      auto am = get_symm<metric_dim>(values, av);
      auto ax = get_vector<mesh_dim>(coords, av);
      auto vec = ax - x;
      auto metric_dist = metric_length(am, vec);
      auto factor = metric_eigenvalue_from_length(1.0 + metric_dist * max_rate);
      auto limiter = am * factor;
      auto limited = intersect_metrics(m, limiter);
      m = limited;
    }
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f, "limit_metric_gradation");
  values = Reals(out);
  values = mesh->sync_array(VERT, values, symm_ncomps(metric_dim));
  return values;
}

static Reals limit_gradation_once(Mesh* mesh, Reals values, Real max_rate) {
  auto metric_dim = get_metrics_dim(mesh->nverts(), values);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return limit_gradation_once_tmpl<3, 3>(mesh, values, max_rate);
  } else if (mesh->dim() == 2 && metric_dim == 2) {
    return limit_gradation_once_tmpl<2, 2>(mesh, values, max_rate);
  } else if (mesh->dim() == 3 && metric_dim == 1) {
    return limit_gradation_once_tmpl<3, 1>(mesh, values, max_rate);
  } else if (mesh->dim() == 2 && metric_dim == 1) {
    return limit_gradation_once_tmpl<2, 1>(mesh, values, max_rate);
  } else if (mesh->dim() == 1) {
    return limit_gradation_once_tmpl<1, 1>(mesh, values, max_rate);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals limit_metric_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol, bool verbose) {
  OMEGA_H_TIME_FUNCTION;
  OMEGA_H_CHECK(mesh->owners_have_all_upward(VERT));
  OMEGA_H_CHECK(max_rate > 0.0);
  auto comm = mesh->comm();
  Reals values2 = values;
  Int i = 0;
  do {
    values = values2;
    values2 = limit_gradation_once(mesh, values, max_rate);
    ++i;
    if (verbose && can_print(mesh) && i > 0 && i % 50 == 0) {
      std::cout << "warning: gradation limiting is up to step " << i << '\n';
    }
  } while (!comm->reduce_and(are_close(values, values2, tol)));
  if (verbose && can_print(mesh)) {
    std::cout << "limited gradation in " << i << " steps\n";
  }
  return values2;
}

template <Int metric_dim>
Reals project_metrics_dim(Mesh* mesh, Reals e2m) {
  auto e_linear = linearize_metrics(mesh->nelems(), e2m);
  auto v2e = mesh->ask_up(VERT, mesh->dim());
  auto v_metrics_w = Write<Real>(mesh->nverts() * symm_ncomps(metric_dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto vm = zero_matrix<metric_dim, metric_dim>();
    Int n = 0;
    for (auto ve = v2e.a2ab[v]; ve < v2e.a2ab[v + 1]; ++ve) {
      auto e = v2e.ab2b[ve];
      auto em = get_symm<metric_dim>(e2m, e);
      average_metric_contrib(vm, n, em, true);
    }
    vm = average_metric_finish(vm, n, true);
    set_symm(v_metrics_w, v, vm);
  };
  parallel_for(mesh->nverts(), f);
  auto v_metrics = Reals(v_metrics_w);
  return mesh->sync_array(VERT, v_metrics, symm_ncomps(metric_dim));
}

Reals project_metrics(Mesh* mesh, Reals e2m) {
  auto metric_dim = get_metrics_dim(mesh->nelems(), e2m);
  if (metric_dim == 3)
    return project_metrics_dim<3>(mesh, e2m);
  else if (metric_dim == 2)
    return project_metrics_dim<2>(mesh, e2m);
  else if (metric_dim == 1)
    return project_metrics_dim<1>(mesh, e2m);
  else
    OMEGA_H_NORETURN(Reals());
}

Reals smooth_metric_once(Mesh* mesh, Reals v2m, bool has_degen) {
  return project_metrics(
      mesh, get_mident_metrics(mesh, mesh->dim(), v2m, has_degen));
}

template <Int dim>
Reals element_implied_length_metrics_dim(Mesh* mesh) {
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto sizes = mesh->ask_sizes();
  auto out = Write<Real>(mesh->nelems() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO e) {
    auto v = gather_verts<dim + 1>(ev2v, e);
    auto p = gather_vectors<dim + 1, dim>(coords, v);
    auto m = element_implied_metric(p);
    set_symm(out, e, m);
  };
  parallel_for(mesh->nelems(), f, "element_implied_length_metrics");
  return out;
}

Reals get_element_implied_length_metrics(Mesh* mesh) {
  if (mesh->dim() == 3) return element_implied_length_metrics_dim<3>(mesh);
  if (mesh->dim() == 2) return element_implied_length_metrics_dim<2>(mesh);
  if (mesh->dim() == 1) return element_implied_length_metrics_dim<1>(mesh);
  OMEGA_H_NORETURN(Reals());
}

Reals get_pure_implied_metrics(Mesh* mesh) {
  return project_metrics(mesh, get_element_implied_length_metrics(mesh));
}

/* These are completely empirical estimates of the volume of
   an element in a "real" unit-edge-length mesh */
static constexpr Real typical_unit_simplex_size(Int dim) {
  return (dim == 3 ? 0.0838934100219 : (dim == 2 ? 0.377645136635 : 1.0));
}

static constexpr Real get_typical_over_perfect_size(Int dim) {
  return typical_unit_simplex_size(dim) / equilateral_simplex_size(dim);
}

static Reals get_element_implied_size_metrics(Mesh* mesh) {
  auto length_metrics = get_element_implied_length_metrics(mesh);
  auto typical_over_perfect = get_typical_over_perfect_size(mesh->dim());
  auto size_scalar = typical_over_perfect;
  auto metric_scalar = power(size_scalar, 2, mesh->dim());
  return multiply_each_by(length_metrics, metric_scalar);
}

Reals get_implied_metrics(Mesh* mesh) {
  begin_code("get_implied_metrics");
  auto out = project_metrics(mesh, get_element_implied_size_metrics(mesh));
  end_code();
  return out;
}

Reals get_pure_implied_isos(Mesh* mesh) {
  auto metrics = get_pure_implied_metrics(mesh);
  return apply_isotropy(mesh->nverts(), metrics, OMEGA_H_ISO_SIZE);
}

Reals get_implied_isos(Mesh* mesh) {
  OMEGA_H_TIME_FUNCTION;
  auto metrics = get_implied_metrics(mesh);
  metrics = apply_isotropy(mesh->nverts(), metrics, OMEGA_H_ISO_SIZE);
  return metrics;
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
static OMEGA_H_INLINE_BIG Tensor<dim> metric_from_hessian(
    Tensor<dim> const hessian, Real const eps) {
  auto const ed = decompose_eigen(hessian);
  auto const r = ed.q;
  auto const l = ed.l;
  constexpr auto c_num = square(dim);
  constexpr auto c_denom = 2 * square(dim + 1);
  decltype(ed.l) tilde_l;
  for (Int i = 0; i < dim; ++i) {
    tilde_l[i] = (c_num * std::abs(l[i])) / (c_denom * eps);
  }
  return compose_eigen(r, tilde_l);
}

template <Int dim>
Reals metric_from_hessians_dim(Reals hessians, Real eps) {
  auto ncomps = symm_ncomps(dim);
  auto n = divide_no_remainder(hessians.size(), ncomps);
  auto out = Write<Real>(n * ncomps);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto hess = get_symm<dim>(hessians, i);
    auto m = metric_from_hessian(hess, eps);
    set_symm(out, i, m);
  };
  parallel_for(n, f, "metric_from_hessians");
  return out;
}

Reals get_hessian_metrics(Int dim, Reals hessians, Real eps) {
  OMEGA_H_CHECK(eps > 0.0);
  if (dim == 3) return metric_from_hessians_dim<3>(hessians, eps);
  if (dim == 2) return metric_from_hessians_dim<2>(hessians, eps);
  if (dim == 1) return metric_from_hessians_dim<1>(hessians, eps);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
static OMEGA_H_INLINE_BIG Tensor<dim> metric_from_gradient(
    Vector<dim> const grad, Real const eps) {
  auto const grad_norm_sq = norm_squared(grad);
  constexpr auto c_num = square(dim);
  constexpr auto c_denom = square(2 * (dim + 1));
  auto const l = (c_num * grad_norm_sq) / (c_denom * square(eps));
  if (l < EPSILON) return zero_matrix<dim, dim>();
  auto const grad_norm = std::sqrt(grad_norm_sq);
  auto const dir = grad / grad_norm;
  return outer_product(dir, dir) * l;
}

template <Int dim>
Reals metric_from_gradients_dim(Reals gradients, Real eps) {
  auto n = divide_no_remainder(gradients.size(), dim);
  auto out = Write<Real>(n * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO const i) {
    auto const grad = get_vector<dim>(gradients, i);
    auto const m = metric_from_gradient(grad, eps);
    set_symm(out, i, m);
  };
  parallel_for(n, f, "metric_from_gradients");
  return out;
}

Reals get_gradient_metrics(Int dim, Reals gradients, Real eps) {
  OMEGA_H_CHECK(eps > 0.0);
  if (dim == 3) return metric_from_gradients_dim<3>(gradients, eps);
  if (dim == 2) return metric_from_gradients_dim<2>(gradients, eps);
  if (dim == 1) return metric_from_gradients_dim<1>(gradients, eps);
  OMEGA_H_NORETURN(Reals());
}

/* this algorithm creates degenerate metrics that only specify size in
   the tangential direction to the sharp curves on the mesh surface */
template <Int dim>
void get_curve_curvature_metrics(
    SurfaceInfo surface_info, Real segment_angle, Write<Real> out) {
  auto f = OMEGA_H_LAMBDA(LO const curv_vert) {
    auto const k = surface_info.curv_vert_curvatures[curv_vert];
    auto const t = get_vector<dim>(surface_info.curv_vert_tangents, curv_vert);
    auto const ew = square(k / segment_angle);
    auto const m = outer_product(t, ew * t);  // t * ew * transpose(t)
    auto const vert = surface_info.curv_vert2vert[curv_vert];
    set_symm(out, vert, m);
  };
  parallel_for(surface_info.curv_vert2vert.size(), std::move(f));
}

Reals get_curvature_metrics(Mesh* mesh, Real segment_angle) {
  auto surface_info = get_surface_info(mesh);
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(mesh->dim()), 0.0);
  if (mesh->dim() == 3) {
    /* this algorithm creates degenerate metrics that only specify size in
       the two tangential directions to mesh surfaces */
    auto f = OMEGA_H_LAMBDA(LO const surf_vert) {
      auto const II = get_symm<2>(surface_info.surf_vert_IIs, surf_vert);
      auto const II_decomp = decompose_eigen(II);
      Vector<2> m_ews;
      for (Int i = 0; i < 2; ++i) {
        m_ews[i] = square(II_decomp.l[i] / segment_angle);
      }
      auto const n = get_vector<3>(surface_info.surf_vert_normals, surf_vert);
      auto const frame = form_ortho_basis(n);
      Matrix<3, 2> surf_frame_t;
      surf_frame_t[0] = frame[1];
      surf_frame_t[1] = frame[2];
      auto const surf_frame = transpose(surf_frame_t);
      auto const m_q_inv = II_decomp.q * surf_frame;
      auto const m_q = pseudo_invert(m_q_inv);
      auto const m = m_q * diagonal(m_ews) * m_q_inv;
      auto const vert = surface_info.surf_vert2vert[surf_vert];
      set_symm(out, vert, m);
    };
    parallel_for(surface_info.surf_vert2vert.size(), std::move(f));
    get_curve_curvature_metrics<3>(surface_info, segment_angle, out);
  } else if (mesh->dim() == 2) {
    get_curve_curvature_metrics<2>(surface_info, segment_angle, out);
  }
  return out;
}

/* The algorithms below are for scaling a size field such that
 * adapting based on that size field will result in a certain specified
 * number of elements.
 *
 * Much of the inspiration came from Section 2.7 of:
 * Pain, C. C., et al.
 * "Tetrahedral mesh optimisation and adaptivity for
 *  steady-state and transient finite element calculations."
 * Computer Methods in Applied Mechanics and Engineering
 * 190.29 (2001): 3771-3796.
 */

template <Int mesh_dim, Int metric_dim>
Reals get_complexity_per_elem_tmpl(Mesh* mesh, Reals v2m) {
  auto const elems2verts = mesh->ask_elem_verts();
  auto const coords = mesh->coords();
  auto const out_w = Write<Real>(mesh->nelems());
  auto const elem_metrics = get_mident_metrics(mesh, mesh->dim(), v2m);
  auto f = OMEGA_H_LAMBDA(LO const e) {
    auto const v = gather_verts<mesh_dim + 1>(elems2verts, e);
    auto const p = gather_vectors<mesh_dim + 1, mesh_dim>(coords, v);
    auto const b = simplex_basis<mesh_dim, mesh_dim>(p);
    auto const real_volume = simplex_size_from_basis(b);
    auto const m = get_symm<metric_dim>(elem_metrics, e);
    auto const sqrt_metric_det =
        power<mesh_dim, 2 * metric_dim>(determinant(m));
    out_w[e] = real_volume * sqrt_metric_det;
  };
  parallel_for(mesh->nelems(), std::move(f));
  return Reals(out_w);
}

Reals get_complexity_per_elem(Mesh* mesh, Reals v2m) {
  if (v2m.size() == 0) return Reals({});
  auto metric_dim = get_metrics_dim(mesh->nverts(), v2m);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return get_complexity_per_elem_tmpl<3, 3>(mesh, v2m);
  } else if (mesh->dim() == 2 && metric_dim == 2) {
    return get_complexity_per_elem_tmpl<2, 2>(mesh, v2m);
  } else if (mesh->dim() == 3 && metric_dim == 1) {
    return get_complexity_per_elem_tmpl<3, 1>(mesh, v2m);
  } else if (mesh->dim() == 2 && metric_dim == 1) {
    return get_complexity_per_elem_tmpl<2, 1>(mesh, v2m);
  } else if (mesh->dim() == 1) {
    return get_complexity_per_elem_tmpl<1, 1>(mesh, v2m);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals get_nelems_per_elem(Mesh* mesh, Reals v2m) {
  auto complexity = get_complexity_per_elem(mesh, v2m);
  return multiply_each_by(
      complexity, 1.0 / typical_unit_simplex_size(mesh->dim()));
}

Real get_complexity(Mesh* mesh, Reals v2m) {
  auto complexity_per_elem = get_complexity_per_elem(mesh, v2m);
  return repro_sum_owned(mesh, mesh->dim(), complexity_per_elem);
}

Real get_expected_nelems_from_complexity(Real complexity, Int dim) {
  return complexity / typical_unit_simplex_size(dim);
}

Real get_expected_nelems(Mesh* mesh, Reals v2m) {
  auto complexity = get_complexity(mesh, v2m);
  return get_expected_nelems_from_complexity(complexity, mesh->dim());
}

Real get_metric_scalar_for_nelems(
    Int elem_dim, Real expected_nelems, Real target_nelems) {
  auto size_scal = target_nelems / expected_nelems;
  auto metric_scal = power(size_scal, 2, elem_dim);
  return metric_scal;
}

Real get_metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems) {
  auto nelems = get_expected_nelems(mesh, v2m);
  auto metric_scal =
      get_metric_scalar_for_nelems(mesh->dim(), nelems, target_nelems);
  return metric_scal;
}

template <Int dim>
Reals intersect_metrics_dim(Reals a, Reals b) {
  auto n = divide_no_remainder(a.size(), symm_ncomps(dim));
  OMEGA_H_CHECK(a.size() == b.size());
  auto c = Write<Real>(n * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto am = get_symm<dim>(a, i);
    auto bm = get_symm<dim>(b, i);
    auto cm = intersect_metrics(am, bm);
    set_symm(c, i, cm);
  };
  parallel_for(n, f, "intersect_metrics");
  return c;
}

Reals intersect_metrics(LO nmetrics, Reals a, Reals b) {
  auto dim = get_metrics_dim(nmetrics, a);
  if (dim == 1) return intersect_metrics_dim<1>(a, b);
  if (dim == 2) return intersect_metrics_dim<2>(a, b);
  if (dim == 3) return intersect_metrics_dim<3>(a, b);
  OMEGA_H_NORETURN(Reals());
}

template <Int new_dim>
Reals metrics_from_isos_dim(Reals isos) {
  auto n = isos.size();
  Write<Real> new_symms(n * symm_ncomps(new_dim));
  auto f = OMEGA_H_LAMBDA(Int i) {
    set_symm(new_symms, i, diagonal(fill_vector<new_dim>(isos[i])));
  };
  parallel_for(n, f, "metrics_from_isos");
  return new_symms;
}

Reals metrics_from_isos(Int new_dim, Reals isos) {
  if (new_dim == 1) return isos;
  if (new_dim == 2) return metrics_from_isos_dim<2>(isos);
  if (new_dim == 3) return metrics_from_isos_dim<3>(isos);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals get_size_isos_dim(Reals metrics) {
  auto n = divide_no_remainder(metrics.size(), symm_ncomps(dim));
  auto out = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto m = get_symm<dim>(metrics, i);
    out[i] = root<dim>(determinant(m));
  };
  parallel_for(n, f, "get_size_isos");
  return out;
}

static Reals get_size_isos(Int dim, Reals metrics) {
  switch (dim) {
    case 3:
      return get_size_isos_dim<3>(metrics);
    case 2:
      return get_size_isos_dim<2>(metrics);
    case 1:
      return get_size_isos_dim<1>(metrics);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals apply_isotropy(LO nmetrics, Reals metrics, Omega_h_Isotropy isotropy) {
  if (nmetrics == 0) return metrics;
  auto dim = get_metrics_dim(nmetrics, metrics);
  switch (isotropy) {
    case OMEGA_H_ANISOTROPIC:
      return metrics;
    case OMEGA_H_ISO_LENGTH:
      return get_max_eigenvalues(dim, metrics);
    case OMEGA_H_ISO_SIZE:
      return get_size_isos(dim, metrics);
  }
  OMEGA_H_NORETURN(Reals());
}

Reals isos_from_lengths(Reals h) {
  auto out = Write<Real>(h.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    out[i] = metric_eigenvalue_from_length(h[i]);
  };
  parallel_for(h.size(), f);
  return out;
}

Reals lengths_from_isos(Reals l) {
  auto out = Write<Real>(l.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    out[i] = metric_length_from_eigenvalue(l[i]);
  };
  parallel_for(l.size(), f);
  return out;
}

}  // end namespace Omega_h
