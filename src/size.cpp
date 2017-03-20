#include "size.hpp"

#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_map.hpp"
#include "eigen.hpp"
#include "project.hpp"
#include "quality.hpp"
#include "surface.hpp"

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
  auto metrics = mesh->get_array(VERT, "metric");
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  if (mesh->dim() == 3 && metric_dim == 3) {
    return measure_edges_tmpl<MetricEdgeLengths<3, 3>>(mesh, a2e);
  }
  if (mesh->dim() == 2 && metric_dim == 2) {
    return measure_edges_tmpl<MetricEdgeLengths<2, 2>>(mesh, a2e);
  }
  if (mesh->dim() == 3 && metric_dim == 1) {
    return measure_edges_tmpl<MetricEdgeLengths<3, 1>>(mesh, a2e);
  }
  if (mesh->dim() == 2 && metric_dim == 1) {
    return measure_edges_tmpl<MetricEdgeLengths<2, 1>>(mesh, a2e);
  }
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
  return project_metrics(mesh, element_implied_metrics(mesh));
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
  auto msl_obj = MeanSquaredLength(mesh, v2m);
  auto ev2v = mesh->ask_elem_verts();
  auto coords = mesh->coords();
  auto out_w = Write<Real>(mesh->nelems());
  auto f = LAMBDA(LO e) {
    auto eev2v = gather_verts<mesh_dim + 1>(ev2v, e);
    auto eev2x = gather_vectors<mesh_dim + 1, mesh_dim>(coords, eev2v);
    auto basis = simplex_basis<mesh_dim, mesh_dim>(eev2x);
    auto edge_vectors = element_edge_vectors(eev2x, basis);
    auto msl = msl_obj.template get<metric_dim>(e, edge_vectors);
    out_w[e] = power<dim, 2>(msl);
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

Reals get_curvature_isos(Mesh* mesh, Real segment_angle) {
  auto vert_curvatures = get_vert_curvatures(mesh);
  auto max_radius = max_size / segment_angle;
  auto out = Write<Real>(mesh->nverts());
  auto f = LAMBDA(LO v) {
    auto curvature = vert_curvatures[v];
    auto l = square(curvature / segment_angle);
    out[v] = l;
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

Reals get_gradient_isos(
    Mesh* mesh, Real error_bound, Reals scalar_field) {
  auto gradients = recover_gradients(mesh, scalar_field);
  auto norms = get_vector_norms(gradients, mesh->dim());
  auto out = Write<Real>(mesh->nverts());
  auto f = LAMBDA(LO v) {
    auto u_dot = norms[v];
    auto l = square(u_dot / error_bound);
    out[v] = l;
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

}  // end namespace Omega_h
