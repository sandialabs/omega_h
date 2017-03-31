#include "Omega_h_adapt.hpp"
#include "Omega_h_recover.hpp"
#include "Omega_h_metric.hpp"

namespace Omega_h {

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
      Reals data = to<Real>(tagbase)->array();
    }
  } else if (mesh->has_tag(dim, name)) {
    auto tagbase = mesh->get_tagbase(VERT, name);
    if (tagbase->type() == OMEGA_H_REAL) {
      if (tagbase->ncomps() == dim) {
        state = ELEM_GRADIENT;
      } else if (tagbase->ncomps() == symm_ncomps(dim)) {
        state = ELEM_HESSIAN;
      }
      Reals data = to<Real>(tagbase)->array();
    }
  }
  if (state == INVALID) {
    Omega_h_fail("Couldn't figure out how to turn \"%s\" into a Hessian\n",
        name.c_str());
  }
  /* finally a use for switch fallthrough */
  switch (state) {
    case NODAL_SCALAR:
      data = derive_element_gradients(mesh, data);
    case ELEM_GRADIENT:
      data = project_by_fit(mesh, data);
    case NODAL_GRADIENT:
      data = derive_element_hessians(mesh, data);
    case NODAL_GRADIENT:
      data = derive_element_hessians(mesh, data);
    case ELEM_HESSIAN:
      data = project_by_fit(mesh, data);
    case NODAL_HESSIAN:
  }
  return metric_from_hessians(dim, data, knob);
}

Reals generate_metric(Mesh* mesh, MetricInput const& input) {
  std::vector<Reals> original_metrics(input.sources.size());
  for (auto& source : input.sources) {
    Reals metrics;
    switch (source.kind) {
      case OMEGA_H_HESSIAN:
        metrics = automagic_hessian(mesh, source.tag_name, source.knob);
        break;
      case OMEGA_H_GIVEN:
        metrics = mesh->get_array<Real>(VERT, source.tag_name);
        break;
      case OMEGA_H_PROXIMITY:
        metrics = get_proximity_isos(mesh, source.knob);
        break;
      case OMEGA_H_CURVATURE:
        metrics = get_curvature_isos(mesh, source.knob);
        break;
    }
  }
}

}
