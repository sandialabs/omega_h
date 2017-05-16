#ifndef OMEGA_H_MIGRATE_HPP
#define OMEGA_H_MIGRATE_HPP

#include <Omega_h_adj.hpp>
#include <Omega_h_dist.hpp>

namespace Omega_h {

class Mesh;

/* create arrays mapping uses of (low_dim) entities by
   (high_dim) entities to their (low_dim) owners */
Remotes form_down_use_owners(Mesh* mesh, Int high_dim, Int low_dim);

/* Given a parallel map (Dist) from
   (new local copies of uses of low entities by high entities)
   to (old owner copies of low entities), this function finds
   use copies that reside on the same MPI rank and have the same
   (old owner copy of low entity), ignoring what their (high entity) is.
   It then removes such duplicates, and what remains is
   a set of uses that are uniquely defined by (old low owner)
   and MPI rank, meaning they are exactly the set of new (low entity)
   copies.
   It then returns a parallel map (Dist) from the (new low copies),
   i.e. the unique uses, to the (old low owner copies).
   This function is the one responsible for the ordering
   of new entities (that are not elements).
   It will index them in order of ascending global number. */
Dist get_new_copies2old_owners(Dist uses2old_owners, GOs old_owner_globals);

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
    Dist new_ents2old_owners, Dist old_owners2new_ents, Omega_h_Parting mode);

void migrate_mesh(
    Mesh* mesh, Dist new_elems2old_owners, Omega_h_Parting mode, bool verbose);

}  // end namespace Omega_h

#endif
