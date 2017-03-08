#ifndef OMEGA_H_SIZE_FIELD_HPP
#define OMEGA_H_SIZE_FIELD_HPP

#include <string>

#include <Omega_h_adapt.hpp>
#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

Reals find_implied_size(Mesh* mesh);
Reals find_implied_metric(Mesh* mesh);
void axes_from_metric_field(
    Mesh* mesh, std::string const& metric_name, std::string const& axis_prefix);
Reals limit_size_field_gradation(
    Mesh* mesh, Reals values, Real max_rate, Real tol = 1e-3);
Reals expected_elems_per_elem_iso(Mesh* mesh, Reals v2h);
Reals expected_elems_per_elem_metric(Mesh* mesh, Reals v2m);
Real size_scalar_for_nelems(Mesh* mesh, Reals v2h, Real target_nelems);
Real metric_scalar_for_nelems(Mesh* mesh, Reals v2m, Real target_nelems);
Reals smooth_metric_once(Mesh* mesh, Reals v2m);
Reals smooth_isos_once(Mesh* mesh, Reals v2h);
Reals get_curvature_isos(Mesh* mesh, Real segment_angle, Real max_size);
Reals get_gradient_isos(
    Mesh* mesh, Real error_bound, Real max_size, Reals scalar_field);
Reals clamp_deforming_isos(Mesh* mesh, Reals isos, Real min_size,
    Real max_interior_size, Real max_boundary_size);

Reals recover_hessians(Mesh* mesh, Reals vert_values);
Reals metric_from_hessians(Int dim, Reals hessians, Real eps, Real hmax);
Reals metric_for_nelems_from_hessians(
    Mesh* mesh, Real target_nelems, Real tolerance, Reals hessians, Real hmax);

}

#endif
