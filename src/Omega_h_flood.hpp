#ifndef OMEGA_H_FLOOD_HPP
#define OMEGA_H_FLOOD_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_adapt.hpp>

namespace Omega_h {

class Mesh;

Bytes mark_floodable_elements(Mesh* mesh);
Bytes mark_flood_seeds(Mesh* mesh, AdaptOpts const& opts,
    Bytes elems_can_flood);
Bytes mark_seeded_flood_zones(Mesh* mesh, Bytes elems_can_flood, Bytes elems_are_seeded);
void flood_element_variables(Mesh* mesh, 
    Bytes elems_should_flood,
    Reals elem_densities,
    Read<I32>* p_elem_flood_class_ids,
    Reals* p_elem_flood_densities);
void flood_class_ids(Mesh* mesh, Int ent_dim);
void flood_classification(Mesh* mesh, Bytes elems_did_flood);

}  // end namespace Omega_h

#endif
