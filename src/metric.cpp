#include "metric.hpp"

#include "access.hpp"
#include "loop.hpp"

namespace Omega_h {

template <Int sdim, Int edim>
static Reals average_metric_tmpl(Mesh* mesh, LOs a2e, Reals v2m) {
  auto na = a2e.size();
  Write<Real> out(na * symm_dofs(sdim));
  auto ev2v = mesh->ask_verts_of(edim);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<edim + 1>(ev2v, e);
    auto ms = gather_symms<edim + 1, sdim>(v2m, v);
    auto m = average_metrics(ms);
    set_symm(out, a, m);
  };
  parallel_for(na, f);
  return out;
}

Reals average_metric(Mesh* mesh, Int ent_dim, LOs entities, Reals v2m) {
  if (mesh->dim() == 3) {
    if (ent_dim == 3) {
      return average_metric_tmpl<3, 3>(mesh, entities, v2m);
    } else if (ent_dim == 2) {
      return average_metric_tmpl<3, 2>(mesh, entities, v2m);
    } else {
      CHECK(ent_dim == 1);
      return average_metric_tmpl<3, 1>(mesh, entities, v2m);
    }
  } else {
    CHECK(mesh->dim() == 2);
    if (ent_dim == 2) {
      return average_metric_tmpl<2, 2>(mesh, entities, v2m);
    } else {
      CHECK(ent_dim == 1);
      return average_metric_tmpl<2, 1>(mesh, entities, v2m);
    }
  }
}

template <Int dim>
Reals interpolate_metrics(Reals a, Reals b, Real t) {
  CHECK(a.size() == b.size());
  CHECK(a.size() % symm_dofs(dim) == 0);
  auto n = a.size() / symm_dofs(dim);
  auto out = Write<Real>(n * symm_dofs(dim));
  auto f = LAMBDA(LO i) {
    auto am = get_symm<dim>(a, i);
    auto bm = get_symm<dim>(b, i);
    auto cm = interpolate_metrics(am, bm, t);
    set_symm(out, i, cm);
  };
  parallel_for(n, f);
  return out;
}

Reals interpolate_metrics(Int dim, Reals a, Reals b, Real t) {
  if (dim == 3) return interpolate_metrics<3>(a, b, t);
  CHECK(dim == 2);
  return interpolate_metrics<2>(a, b, t);
}

template <Int dim>
Reals linearize_metrics_dim(Reals metrics) {
  auto n = metrics.size() / symm_dofs(dim);
  auto out = Write<Real>(n * dim * dim);
  auto f = LAMBDA(LO i) {
    set_matrix(out, i, linearize_metric(get_symm<dim>(metrics, i)));
  };
  parallel_for(n, f);
  return out;
}

template <Int dim>
Reals delinearize_metrics_dim(Reals lms) {
  auto n = lms.size() / square(dim);
  auto out = Write<Real>(n * symm_dofs(dim));
  auto f = LAMBDA(LO i) {
    set_symm(out, i, delinearize_metric(get_matrix<dim>(lms, i)));
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
  CHECK(linear_metrics.size() % square(dim) == 0);
  if (dim == 3) return delinearize_metrics_dim<3>(linear_metrics);
  if (dim == 2) return delinearize_metrics_dim<2>(linear_metrics);
  NORETURN(Reals());
}

}  // end namespace Omega_h
