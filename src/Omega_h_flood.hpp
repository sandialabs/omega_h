#ifndef OMEGA_H_FLOOD_HPP
#define OMEGA_H_FLOOD_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

class Mesh;

struct FloodOpts {
  std::string density_name;
  Real min_length;
  Real min_angle;
};

Bytes mark_flood_zones(Mesh* mesh, FloodOpts const& opts);
void flood_element_variables(Mesh* mesh, Bytes elems_should_flood,
    Reals elem_densities, Read<I32>* p_elem_flood_class_ids,
    Reals* p_elem_flood_densities);
void flood_class_ids(Mesh* mesh, Int ent_dim);
void flood_classification(Mesh* mesh, Bytes elems_did_flood);
void flood(Mesh* mesh, FloodOpts const& opts);

}  // end namespace Omega_h

#endif
