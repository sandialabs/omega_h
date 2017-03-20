#ifndef OMEGA_H_SIZE_FIELD_HPP
#define OMEGA_H_SIZE_FIELD_HPP

#include <string>

#include <Omega_h_adapt.hpp>
#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

Reals find_implied_isos(Mesh* mesh);
Reals find_implied_metric(Mesh* mesh);
void axes_from_metric_field(
    Mesh* mesh, std::string const& metric_name, std::string const& axis_prefix);
Reals limit_size_field_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol = 1e-3);
Reals expected_elems_per_elem(Mesh* mesh, Reals v2m);
Real metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems);
Reals smooth_metric_once(Mesh* mesh, Reals v2m);
Reals smooth_isos_once(Mesh* mesh, Reals v2h);
Reals get_curvature_isos(Mesh* mesh, Real segment_angle);
Reals get_gradient_isos(
    Mesh* mesh, Real error_bound, Reals scalar_field);
Reals get_aniso_zz_metric(
    Mesh* mesh, Reals elem_gradients, Real error_bound, Real max_size);

Reals recover_hessians(Mesh* mesh, Reals vert_values);
Reals metric_from_hessians(Int dim, Reals hessians, Real eps, Real hmax);
Reals metric_for_nelems_from_hessians(
    Mesh* mesh, Real target_nelems, Real tolerance, Reals hessians, Real hmax);

Reals derive_element_gradients(Mesh* mesh, Reals vert_values);
Reals derive_element_hessians(Mesh* mesh, Reals vert_gradients);
Reals recover_gradients(Mesh* mesh, Reals vert_values);
Reals recover_hessians_from_gradients(Mesh* mesh, Reals vert_gradients);

Reals project_metrics(Mesh* mesh, Reals e2m);

}  // namespace Omega_h

#endif
