#ifndef TRANSFER_HPP
#define TRANSFER_HPP

#include "internal.hpp"

namespace osh {

void transfer_refine(Mesh* old_mesh, Mesh* new_mesh, LOs keys2edges,
                     LOs keys2midverts, Int prod_dim, LOs keys2prods,
                     LOs prods2new_ents, LOs same_ents2old_ents,
                     LOs same_ents2new_ents);

void transfer_coarsen(Mesh* old_mesh, Mesh* new_mesh, LOs keys2verts,
                      Adj keys2doms, Int prod_dim, LOs prods2new_ents,
                      LOs same_ents2old_ents, LOs same_ents2new_ents);

void transfer_swap(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim, LOs keys2edges,
                   LOs keys2prods, LOs prods2new_ents, LOs same_ents2old_ents,
                   LOs same_ents2new_ents);

void transfer_copy(Mesh* old_mesh, Mesh* new_mesh, Int prod_dim);

}  // end namespace osh

#endif
