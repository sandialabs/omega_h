#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_profile.hpp"
#include "Omega_h_recover.hpp"
#include "Omega_h_timer.hpp"

#include <iostream>

namespace Omega_h {

MetricSource::MetricSource(Omega_h_Source type_, Real knob_,
    std::string const& tag_name_, Omega_h_Isotropy isotropy_,
    Omega_h_Scales scales_)
    : type(type_),
      knob(knob_),
      tag_name(tag_name_),
      isotropy(isotropy_),
      scales(scales_) {}

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
  min_element_count = 1.0;
  element_count_over_relaxation = 1.1;
  nsmoothing_steps = 0;
}

void MetricInput::add_source(MetricSource const& src) {
  sources.push_back(src);
}

static Reals get_variation_metrics(
    Mesh* mesh, Real knob, Int dim, Int ncomps, Reals data) {
  OMEGA_H_CHECK(data.size() == mesh->nents(dim) * ncomps);
  if (ncomps == 1) {
    if (dim == VERT) {
      auto hessians = recover_hessians(mesh, data);
      return get_hessian_metrics(mesh->dim(), hessians, knob);
    } else if (dim == mesh->dim()) {
      auto vert_data = project_by_fit(mesh, data);
      auto vert_grads = recover_gradients(mesh, vert_data);
      return get_gradient_metrics(mesh->dim(), vert_grads, knob);
    }
  } else {
    Reals metrics;
    for (Int comp = 0; comp < ncomps; ++comp) {
      auto comp_data = get_component(data, ncomps, comp);
      auto comp_metrics = get_variation_metrics(mesh, knob, dim, 1, comp_data);
      if (comp) {
        metrics = intersect_metrics(mesh->nverts(), metrics, comp_metrics);
      } else {
        metrics = comp_metrics;
      }
    }
    return metrics;
  }
  OMEGA_H_NORETURN(Reals());
}

Reals get_variation_metrics(Mesh* mesh, std::string const& name, Real knob) {
  Int dim = -1;
  if (mesh->has_tag(VERT, name)) {
    dim = VERT;
  } else if (mesh->has_tag(mesh->dim(), name)) {
    dim = mesh->dim();
  } else {
    Omega_h_fail("tag \"%s\" neither on elements nor nodes\n", name.c_str());
  }
  auto tag = mesh->get_tag<Real>(dim, name);
  auto ncomps = tag->ncomps();
  auto data = tag->array();
  return get_variation_metrics(mesh, knob, dim, ncomps, data);
}

Reals get_derivative_metrics(Mesh* mesh, std::string const& name, Real knob) {
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
  return get_hessian_metrics(dim, data, knob);
}

Reals generate_metrics(Mesh* mesh, MetricInput const& input) {
  begin_code("generate_metrics");
  auto t0 = now();
  if (input.should_limit_lengths) {
    OMEGA_H_CHECK(input.min_length <= input.max_length);
  }
  if (input.should_limit_element_count) {
    OMEGA_H_CHECK(input.min_length <= input.max_length);
  }
  auto n = mesh->nverts();
  if (!input.sources.size()) {
    Omega_h_fail("generate_metrics: no sources given!\n");
  }
  std::vector<Reals> original_metrics;
  Int metric_dim = 1;
  for (auto& source : input.sources) {
    Reals metrics;
    switch (source.type) {
      case OMEGA_H_CONSTANT:
        metrics =
            Reals(mesh->nverts(), metric_eigenvalue_from_length(source.knob));
        break;
      case OMEGA_H_VARIATION:
        metrics = get_variation_metrics(mesh, source.tag_name, source.knob);
        break;
      case OMEGA_H_DERIVATIVE:
        metrics = get_derivative_metrics(mesh, source.tag_name, source.knob);
        break;
      case OMEGA_H_GIVEN:
        metrics = mesh->get_array<Real>(VERT, source.tag_name);
        if (source.knob != 1.0) {
          metrics = multiply_each_by(
              metrics, metric_eigenvalue_from_length(source.knob));
        }
        break;
      case OMEGA_H_IMPLIED:
        metrics = get_implied_metrics(mesh);
        if (source.knob != 1.0) {
          metrics = multiply_each_by(
              metrics, metric_eigenvalue_from_length(source.knob));
        }
        break;
      case OMEGA_H_CURVATURE:
        metrics = get_curvature_metrics(mesh, source.knob);
        break;
    }
    metrics = apply_isotropy(n, metrics, source.isotropy);
    auto source_dim = get_metrics_dim(n, metrics);
    if (mesh->dim() > 1 && source_dim == mesh->dim()) {
      metric_dim = mesh->dim();
    }
    original_metrics.push_back(metrics);
  }
  Real scalar = 1.0;
  Reals metrics;
  Int niters;
  for (niters = 0; true; ++niters) {
    if (niters == 100) {
      Omega_h_fail("Too many element count limiting iterations\n");
    }
    metrics = Reals();
    for (size_t i = 0; i < input.sources.size(); ++i) {
      auto in_metrics = original_metrics[i];
      if (get_metrics_dim(n, in_metrics) == 1) {
        in_metrics = metrics_from_isos(metric_dim, in_metrics);
      }
      if (input.sources[i].scales == OMEGA_H_SCALES) {
        in_metrics = multiply_each_by(in_metrics, scalar);
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
    for (Int i = 0; i < input.nsmoothing_steps; ++i) {
      metrics = smooth_metric_once(mesh, metrics);
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
  }
  auto t1 = now();
  if (input.verbose && mesh->comm()->rank() == 0) {
    if (input.should_limit_element_count) {
      std::cout << "generated metrics in " << niters << " iterations and "
                << (t1 - t0) << " seconds\n";
    } else {
      std::cout << "generated metrics in " << (t1 - t0) << " seconds\n";
    }
  }
  end_code();
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

void add_implied_isos_tag(Mesh* mesh) {
  OMEGA_H_TIME_FUNCTION;
  auto metrics = get_implied_isos(mesh);
  add_metric_tag(mesh, metrics, "metric");
}

void add_implied_metric_based_on_target(Mesh* mesh) {
  auto target_tagbase = mesh->get_tagbase(VERT, "target_metric");
  if (mesh->dim() > 1 && target_tagbase->ncomps() == 1) {
    add_implied_isos_tag(mesh);
  } else {
    add_implied_metric_tag(mesh);
  }
}

}  // namespace Omega_h
