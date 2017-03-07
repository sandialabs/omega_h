#include "metric.hpp"

#include <iostream>

#include "access.hpp"
#include "Omega_h_array_ops.hpp"
#include "host_few.hpp"
#include "loop.hpp"
#include "project.hpp"
#include "size.hpp"

namespace Omega_h {

template <Int sdim, Int edim>
static Reals mident_metrics_tmpl(Mesh* mesh, LOs a2e, Reals v2m) {
  auto na = a2e.size();
  Write<Real> out(na * symm_dofs(sdim));
  auto ev2v = mesh->ask_verts_of(edim);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<edim + 1>(ev2v, e);
    auto ms = gather_symms<edim + 1, sdim>(v2m, v);
    auto m = average_metric(ms);
    set_symm(out, a, m);
  };
  parallel_for(na, f);
  return out;
}

Reals get_mident_metrics(Mesh* mesh, Int ent_dim, LOs entities, Reals v2m) {
  if (mesh->dim() == 3 && ent_dim == 3) {
    return mident_metrics_tmpl<3, 3>(mesh, entities, v2m);
  }
  if (mesh->dim() == 3 && ent_dim == 1) {
    return mident_metrics_tmpl<3, 1>(mesh, entities, v2m);
  }
  if (mesh->dim() == 2 && ent_dim == 2) {
    return mident_metrics_tmpl<2, 2>(mesh, entities, v2m);
  }
  if (mesh->dim() == 2 && ent_dim == 1) {
    return mident_metrics_tmpl<2, 1>(mesh, entities, v2m);
  }
  NORETURN(Reals());
}

Reals interpolate_between_metrics(Int dim, Reals a, Reals b, Real t) {
  auto log_a = linearize_metrics(dim, a);
  auto log_b = linearize_metrics(dim, b);
  auto log_c = interpolate_between(log_a, log_b, t);
  return delinearize_metrics(dim, log_c);
}

template <Int dim>
Reals linearize_metrics_dim(Reals metrics) {
  auto n = metrics.size() / symm_dofs(dim);
  auto out = Write<Real>(n * symm_dofs(dim));
  auto f = LAMBDA(LO i) {
    set_symm(out, i, linearize_metric(get_symm<dim>(metrics, i)));
  };
  parallel_for(n, f);
  return out;
}

template <Int dim>
Reals delinearize_metrics_dim(Reals lms) {
  auto n = lms.size() / symm_dofs(dim);
  auto out = Write<Real>(n * symm_dofs(dim));
  auto f = LAMBDA(LO i) {
    set_symm(out, i, delinearize_metric(get_symm<dim>(lms, i)));
  };
  parallel_for(n, f);
  return out;
}

Reals linearize_metrics(Int dim, Reals metrics) {
  CHECK(metrics.size() % symm_dofs(dim) == 0);
  if (dim == 3) return linearize_metrics_dim<3>(metrics);
  if (dim == 2) return linearize_metrics_dim<2>(metrics);
  NORETURN(Reals());
}

Reals delinearize_metrics(Int dim, Reals linear_metrics) {
  CHECK(linear_metrics.size() % symm_dofs(dim) == 0);
  if (dim == 3) return delinearize_metrics_dim<3>(linear_metrics);
  if (dim == 2) return delinearize_metrics_dim<2>(linear_metrics);
  NORETURN(Reals());
}

template <Int dim>
static HostFew<Reals, dim> axes_from_metrics_dim(Reals metrics) {
  CHECK(metrics.size() % symm_dofs(dim) == 0);
  auto n = metrics.size() / symm_dofs(dim);
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
    mesh->add_tag(VERT, output_prefix + '_' + std::to_string(i), dim,
        OMEGA_H_DONT_TRANSFER, OMEGA_H_DO_OUTPUT, axes[i]);
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
    Matrix<dim, dim> hessian, Real eps, Real hmax) {
  auto ed = decompose_eigen(hessian);
  auto r = ed.q;
  auto l = ed.l;
  constexpr auto c_num = square(dim);
  constexpr auto c_denom = 2 * square(dim + 1);
  decltype(l) tilde_l;
  for (Int i = 0; i < dim; ++i) {
    auto val = (c_num * fabs(l[i])) / (c_denom * eps);
    tilde_l[i] = max2(val, 1. / square(hmax));
  }
  return compose_eigen(r, tilde_l);
}

template <Int dim>
static Reals metric_from_hessians_dim(Reals hessians, Real eps, Real hmax) {
  auto ncomps = symm_dofs(dim);
  CHECK(hessians.size() % ncomps == 0);
  auto n = hessians.size() / ncomps;
  auto out = Write<Real>(n * ncomps);
  auto f = LAMBDA(LO i) {
    auto hess = get_symm<dim>(hessians, i);
    auto m = metric_from_hessian(hess, eps, hmax);
    set_symm(out, i, m);
  };
  parallel_for(n, f);
  return out;
}

Reals metric_from_hessians(Int dim, Reals hessians, Real eps, Real hmax) {
  CHECK(hmax > 0.0);
  CHECK(eps > 0.0);
  if (dim == 3) return metric_from_hessians_dim<3>(hessians, eps, hmax);
  if (dim == 2) return metric_from_hessians_dim<2>(hessians, eps, hmax);
  NORETURN(Reals());
}

Reals metric_for_nelems_from_hessians(
    Mesh* mesh, Real target_nelems, Real tolerance, Reals hessians, Real hmax) {
  CHECK(tolerance > 0);
  CHECK(target_nelems > 0);
  auto dim = mesh->dim();
  Real scalar;
  Reals metric;
  Real eps = 1.0;
  Int niters = 0;
  do {
    metric = metric_from_hessians(dim, hessians, eps, hmax);
    scalar = metric_scalar_for_nelems(mesh, metric, target_nelems);
    eps /= scalar;
    ++niters;
  } while (fabs(scalar - 1.0) > tolerance);
  if (can_print(mesh)) {
    std::cout << "after " << niters << " iterations,"
              << " metric targets " << target_nelems << "*" << scalar
              << " elements\n";
  }
  return metric;
}

/* gradation limiting code: */

template <Int dim>
class IsoGradation {
 public:
  using Value = Real;
  static INLINE Value form_limiter(Value h, Vector<dim> v, Real rate) {
    auto real_dist = norm(v);
    auto metric_dist = real_dist / h;
    return h * (1.0 + metric_dist * rate);
  }
  static INLINE Value intersect(Value a, Value b) { return min2(a, b); }
  static DEVICE Value get(Reals const& a, LO i) { return a[i]; }
  static DEVICE void set(Write<Real> const& a, LO i, Value v) { a[i] = v; }
  enum { ndofs = 1 };
};

template <Int dim>
class AnisoGradation {
 public:
  using Value = Matrix<dim, dim>;
  static INLINE Value form_limiter(Value m, Vector<dim> v, Real rate) {
    auto metric_dist = metric_length(m, v);
    auto decomp = decompose_metric(m);
    decomp.l = decomp.l * (1.0 + metric_dist * rate);
    return compose_metric(decomp.q, decomp.l);
  }
  static INLINE Value intersect(Value a, Value b) {
    return intersect_metrics(a, b);
  }
  static DEVICE Value get(Reals const& a, LO i) { return get_symm<dim>(a, i); }
  static DEVICE void set(Write<Real> const& a, LO i, Value v) {
    set_symm(a, i, v);
  }
  enum { ndofs = symm_dofs(dim) };
};

template <Int dim, template <Int> class Gradation>
static INLINE typename Gradation<dim>::Value limit_size_value_by_adj(
    typename Gradation<dim>::Value m, Vector<dim> x,
    typename Gradation<dim>::Value am, Vector<dim> ax, Real rate) {
  auto limiter = Gradation<dim>::form_limiter(am, ax - x, rate);
  auto limited = Gradation<dim>::intersect(m, limiter);
  return limited;
}

template <Int dim, template <Int> class Gradation>
static Reals limit_size_field_once_by_adj_tmpl(
    Mesh* mesh, Reals values, Real max_rate) {
  using G = Gradation<dim>;
  auto v2v = mesh->ask_star(VERT);
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * G::ndofs);
  auto f = LAMBDA(LO v) {
    auto m = G::get(values, v);
    auto x = get_vector<dim>(coords, v);
    for (auto vv = v2v.a2ab[v]; vv < v2v.a2ab[v + 1]; ++vv) {
      auto av = v2v.ab2b[vv];
      auto am = G::get(values, av);
      auto ax = get_vector<dim>(coords, av);
      m = limit_size_value_by_adj<dim, Gradation>(m, x, am, ax, max_rate);
    }
    G::set(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  values = Reals(out);
  values = mesh->sync_array(VERT, values, G::ndofs);
  return values;
}

static Reals limit_size_field_once_by_adj(
    Mesh* mesh, Reals values, Real max_rate) {
  if (mesh->dim() == 3) {
    if (values.size() == symm_dofs(3) * mesh->nverts()) {
      return limit_size_field_once_by_adj_tmpl<3, AnisoGradation>(
          mesh, values, max_rate);
    }
    if (values.size() == mesh->nverts()) {
      return limit_size_field_once_by_adj_tmpl<3, IsoGradation>(
          mesh, values, max_rate);
    }
  }
  if (mesh->dim() == 2) {
    if (values.size() == symm_dofs(2) * mesh->nverts()) {
      return limit_size_field_once_by_adj_tmpl<2, AnisoGradation>(
          mesh, values, max_rate);
    }
    if (values.size() == mesh->nverts()) {
      return limit_size_field_once_by_adj_tmpl<2, IsoGradation>(
          mesh, values, max_rate);
    }
  }
  NORETURN(Reals());
}

Reals limit_size_field_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol) {
  CHECK(mesh->owners_have_all_upward(VERT));
  CHECK(max_rate > 0.0);
  auto comm = mesh->comm();
  Reals values2 = values;
  Int i = 0;
  do {
    values = values2;
    values2 = limit_size_field_once_by_adj(mesh, values, max_rate);
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
  auto e_linear = linearize_metrics(mesh->dim(), e2m);
  auto v_linear = project_by_average(mesh, e_linear);
  return delinearize_metrics(mesh->dim(), v_linear);
}

Reals smooth_metric_once(Mesh* mesh, Reals v2m) {
  auto e2e = LOs(mesh->nelems(), 0, 1);
  return project_metrics(mesh, get_mident_metrics(mesh, mesh->dim(), e2e, v2m));
}

}  // end namespace Omega_h
