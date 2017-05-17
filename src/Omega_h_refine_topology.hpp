#ifndef OMEGA_H_REFINE_TOPOLOGY_HPP
#define OMEGA_H_REFINE_TOPOLOGY_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_scalar.hpp>

namespace Omega_h {

class Mesh;

void refine_domains_to_pairs(Mesh* mesh, Int dim, LOs keys2edges,
    LOs keys2midverts, LOs old_verts2new_verts, LOs& keys2pairs,
    LOs& pair_verts2verts);

void refine_domains_to_cuts(Mesh* mesh, Int dim, LOs keys2edges,
    LOs keys2midverts, LOs old_verts2new_verts, LOs& keys2cuts,
    LOs& cut_verts2verts);

void combine_pairs_and_cuts(Int ent_dim, LOs keys2cuts, LOs keys2pairs,
    LOs cut_verts2verts, LOs pair_verts2verts, LOs& keys2prods,
    LOs& prod_verts2verts);

void refine_products(Mesh* mesh, Int ent_dim, LOs keys2edges, LOs keys2midverts,
    LOs old_verts2new_verts, LOs& keys2prods, LOs& prod_verts2verts);

/* as it happens, the triangle-of-tet-to-vertices template
   specifies all triangles curling outward
   (such that implicit face derivation orients all surface
    faces outward), while the tet-to-vertices convention
   has the bottom face curling inward for consistency
   with Gmsh, PUMI, Simmetrix, etc.
   rather than break either of those two compatibilities,
   we will remember to flip the triangle here to get
   the right orientation for new tets */
template <Int dim>
struct FlipNewElem;
template <>
struct FlipNewElem<2> {
  template <typename T>
  OMEGA_H_INLINE static void flip(T ev[]) {
    (void)ev;
  }
};
template <>
struct FlipNewElem<3> {
  template <typename T>
  OMEGA_H_INLINE static void flip(T ev[]) {
    swap2(ev[1], ev[2]);
  }
};
template <Int dim, typename T>
OMEGA_H_INLINE void flip_new_elem(T ev[]) {
  FlipNewElem<dim>::flip(ev);
}
template <typename T>
OMEGA_H_INLINE void flip_new_elem(Int dim, T ev[]) {
  if (dim == 3) swap2(ev[1], ev[2]);
}

}  // end namespace Omega_h

#endif
