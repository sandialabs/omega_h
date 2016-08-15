#ifndef OWNERS_HPP
#define OWNERS_HPP

#include "internal.hpp"

namespace osh {

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

/* uses update_ownership() with a linear partitioning of
   global numbers as the old partitioning.
   if global numbers obey good linear arrangement properties,
   the old partition should be decent and so this is an
   effective fallback if only globals are available */

Remotes owners_from_globals(
    CommPtr comm, Read<GO> globals, Read<I32> own_ranks);

template <typename T>
Read<T> reduce_data_to_owners(
    Read<T> copy_data, Dist copies2owners, Int ncomps);

void globals_from_owners(Mesh* new_mesh, Int ent_dim);

#define INST_DECL(T)                                                           \
  extern template Read<T> reduce_data_to_owners(                               \
      Read<T> copy_data, Dist copies2owners, Int ncomps);
INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL

}  // end namespace osh

#endif
