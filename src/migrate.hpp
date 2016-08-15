#ifndef MIGRATE_HPP
#define MIGRATE_HPP

#include "internal.hpp"

namespace osh {

/* create arrays mapping uses of (low_dim) entities by
   (high_dim) entities to their (low_dim) owners */
Remotes form_down_use_owners(Mesh* mesh, Int high_dim, Int low_dim);

/* given a Dist mapping new entity uses to their old owners,
   filter out duplicate uses of the same old owner by the
   same rank, and create a Dist mapping old owners to
   rank-unique new entity copies. */
Dist find_unique_use_owners(Dist uses2old_owners);

/* given a Dist mapping from new entity copies to old owners,
   and one from old owners to new entity uses,
   form the new connectivity array, which maps each use
   to the local index of its rank-unique copy */
LOs form_new_conn(Dist new_ents2old_owners, Dist old_owners2new_uses);

/* given a Dist mapping old owners to new copies of (ent_dim)
   entities, project this to (low_dim) entities in the form
   of the new (ent_dim -> low_dim) adjacency arrays as well
   as the Dist mapping old (low_dim) owners to new copies */
void push_down(Mesh* old_mesh, Int ent_dim, Int low_dim,
    Dist old_owners2new_ents, Adj& new_ents2new_lows,
    Dist& old_low_owners2new_lows);

void push_tags(Mesh const* old_mesh, Mesh* new_mesh, Int ent_dim,
    Dist old_owners2new_ents);

void push_ents(Mesh* old_mesh, Mesh* new_mesh, Int ent_dim,
    Dist new_ents2old_owners, Dist old_owners2new_ents, osh_parting mode);

void migrate_mesh(Mesh* old_mesh, Mesh* new_mesh, Dist new_elems2old_owners,
    osh_parting mode, bool verbose);
void migrate_mesh(Mesh* mesh, Dist new_elems2old_owners, bool verbose);
void migrate_mesh(Mesh* mesh, Remotes new_elems2old_owners, bool verbose);

}  // end namespace osh

#endif
