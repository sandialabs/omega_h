#ifndef OMEGA_H_OWNERS_HPP
#define OMEGA_H_OWNERS_HPP

#include <Omega_h_dist.hpp>
#include <Omega_h_remotes.hpp>

namespace Omega_h {

class Mesh;

/* compute owners for copies of a new partitioning,
   based on a mapping (Dist) from new copies to
   old owners (their owners in the old partitioning).

   each of the old owners is made responsible for all
   its new copies and selecting an owner among them.
   thus, if both the old and new partitioning are "good",
   this should be a highly scalable way to update
   parallel connectivity

   the (own_ranks) argument is optional. It may be
   left uninitialized (Read<I32>(), !own_ranks.exists()),
   in which case the owner rank will be chosen with
   a preference for ranks that have less copies,
   and in the case two ranks have the same number of copies,
   the smallest rank will be chosen.
   if (own_ranks) is specified, they will dictate the ownership
   of the new copies and are expected to be consistent.
 */

Remotes update_ownership(Dist new_ents2old_owners, Read<I32> own_ranks);

template <typename T>
Read<T> reduce_data_to_owners(
    Read<T> copy_data, Dist copies2owners, Int ncomps);

/* computes new global numbers of entities of dimension (ent_dim).
   owned entities on the same MPI rank will be numbered consecutively,
   and all entities owned by MPI rank (i) will be numbered before
   those owned by MPI rank (i + 1). */
GOs globals_from_owners(Mesh* mesh, Int ent_dim);

#define OMEGA_H_INST_DECL(T)                                                   \
  extern template Read<T> reduce_data_to_owners(                               \
      Read<T> copy_data, Dist copies2owners, Int ncomps);
OMEGA_H_INST_DECL(I8)
OMEGA_H_INST_DECL(I32)
OMEGA_H_INST_DECL(I64)
OMEGA_H_INST_DECL(Real)
#undef OMEGA_H_INST_DECL

}  // end namespace Omega_h

#endif
