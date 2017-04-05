#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_control.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_recover.hpp"
#include "Omega_h_timer.hpp"

#include <iostream>

namespace Omega_h {

MetricInput::MetricInput() {
  verbose = true;
  should_limit_lengths = false;
  max_length = 0.0;
  min_length = 0.0;
  should_limit_gradation = false;
  max_gradation_rate = 1.0;
  gradation_convergence_tolerance = 1e-3;
  should_limit_element_count = false;
  max_element_count = 1e6;
  min_element_count = 0.0;
  element_count_over_relaxation = 1.1;
}

Reals automagic_hessian(Mesh* mesh, std::string const& name, Real knob) {
  enum {
    INVALID,
    NODAL_SCALAR,
    ELEM_GRADIENT,
    NODAL_GRADIENT,
    ELEM_HESSIAN,
    NODAL_HESSIAN,
  } state = INVALID;
  auto dim = mesh->dim();
  Reals data;
  if (mesh->has_tag(VERT, name)) {
    auto tagbase = mesh->get_tagbase(VERT, name);
    if (tagbase->type() == OMEGA_H_REAL) {
      if (tagbase->ncomps() == 1) {
        state = NODAL_SCALAR;
      } else if (tagbase->ncomps() == dim) {
        state = NODAL_GRADIENT;
      } else if (tagbase->ncomps() == symm_ncomps(dim)) {
        state = NODAL_HESSIAN;
      }
      data = as<Real>(tagbase)->array();
    }
  } else if (mesh->has_tag(dim, name)) {
    auto tagbase = mesh->get_tagbase(VERT, name);
    if (tagbase->type() == OMEGA_H_REAL) {
      if (tagbase->ncomps() == dim) {
        state = ELEM_GRADIENT;
      } else if (tagbase->ncomps() == symm_ncomps(dim)) {
        state = ELEM_HESSIAN;
      }
      data = as<Real>(tagbase)->array();
    }
  }
  /* finally a use for switch fallthrough */
  switch (state) {
    case INVALID:
      Omega_h_fail("Couldn't figure out how to turn \"%s\" into a Hessian\n",
          name.c_str());
    case NODAL_SCALAR:
      data = derive_element_gradients(mesh, data);
      OMEGA_H_FALLTHROUGH;
    case ELEM_GRADIENT:
      data = project_by_fit(mesh, data);
      OMEGA_H_FALLTHROUGH;
    case NODAL_GRADIENT:
      data = derive_element_hessians(mesh, data);
      OMEGA_H_FALLTHROUGH;
    case ELEM_HESSIAN:
      data = project_by_fit(mesh, data);
      OMEGA_H_FALLTHROUGH;
    case NODAL_HESSIAN:;
  }
  return metric_from_hessians(dim, data, knob);
}

Reals generate_metrics(Mesh* mesh, MetricInput const& input) {
  auto t0 = now();
  if (input.should_limit_lengths) {
    OMEGA_H_CHECK(input.min_length <= input.max_length);
  }
  if (input.should_limit_element_count) {
    OMEGA_H_CHECK(input.min_length <= input.max_length);
  }
  auto n = mesh->nverts();
  if (!input.sources.size()) {
    if (input.should_limit_lengths) {
      return Reals(n, input.max_length);
    } else {
      Omega_h_fail("generate_metric: no sources or limits given!\n");
    }
  }
  std::vector<Reals> original_metrics;
  Int metric_dim = 1;
  for (auto& source : input.sources) {
    Reals metrics;
    switch (source.type) {
      case OMEGA_H_HESSIAN:
        metrics = automagic_hessian(mesh, source.tag_name, source.knob);
        break;
      case OMEGA_H_GIVEN:
        metrics = mesh->get_array<Real>(VERT, source.tag_name);
        break;
      case OMEGA_H_IMPLIED:
        metrics = get_implied_metrics(mesh);
        break;
      case OMEGA_H_PROXIMITY:
        metrics = get_proximity_isos(mesh, source.knob);
        break;
      case OMEGA_H_CURVATURE:
        metrics = get_curvature_isos(mesh, source.knob);
        break;
    }
    if ((metric_dim == 1) && (get_metrics_dim(n, metrics) > 1)) {
      metric_dim = mesh->dim();
    }
    original_metrics.push_back(metrics);
  }
  Real scalar = 1.0;
  Reals metrics;
  Int niters = 0;
  while (true) {
    metrics = Reals();
    for (size_t i = 0; i < input.sources.size(); ++i) {
      auto in_metrics = original_metrics[i];
      in_metrics =
          resize_symms(in_metrics, get_metrics_dim(n, in_metrics), metric_dim);
      if (input.sources[i].should_scale) {
        in_metrics = multiply_each_by(scalar, in_metrics);
      }
      if (input.should_limit_lengths) {
        in_metrics =
            clamp_metrics(n, in_metrics, input.min_length, input.max_length);
      }
      if (i) {
        metrics = intersect_metrics(n, metrics, in_metrics);
      } else {
        metrics = in_metrics;
      }
    }
    if (input.should_limit_gradation) {
      metrics = limit_metric_gradation(mesh, metrics, input.max_gradation_rate,
          input.gradation_convergence_tolerance, input.verbose);
    }
    if (!input.should_limit_element_count) {
      break;
    } else {
      auto nelems = get_expected_nelems(mesh, metrics);
      if (nelems > input.max_element_count) {
        scalar *= get_metric_scalar_for_nelems(
            mesh->dim(), nelems, input.max_element_count);
        scalar /= input.element_count_over_relaxation;
      } else if (nelems < input.min_element_count) {
        scalar *= get_metric_scalar_for_nelems(
            mesh->dim(), nelems, input.min_element_count);
        scalar *= input.element_count_over_relaxation;
      } else {
        break;
      }
    }
    ++niters;
  }
  auto t1 = now();
  add_to_global_timer("generating metrics", t1 - t0);
  if (input.verbose) {
    std::cout << "generated metrics in " << niters << " iterations and "
              << (t1 - t0) << " seconds";
  }
  return metrics;
}

void add_metric_tag(Mesh* mesh, Reals metrics, std::string const& name) {
  auto metric_dim = get_metrics_dim(mesh->nverts(), metrics);
  mesh->add_tag(VERT, name, symm_ncomps(metric_dim), metrics);
}

void generate_metric_tag(Mesh* mesh, MetricInput const& input) {
  auto metrics = generate_metrics(mesh, input);
  add_metric_tag(mesh, metrics, "metric");
}

void generate_target_metric_tag(Mesh* mesh, MetricInput const& input) {
  auto metrics = generate_metrics(mesh, input);
  add_metric_tag(mesh, metrics, "target_metric");
}

void add_implied_metric_tag(Mesh* mesh) {
  auto metrics = get_implied_metrics(mesh);
  add_metric_tag(mesh, metrics, "metric");
}

}  // namespace Omega_h
