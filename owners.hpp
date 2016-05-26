/* given the global numbers of local copies
   of a partitioned set of objects, construct
   the mapping from each copy to a unique
   owner copy.

   the (own_ranks) argument may be an uninitialized
   array (!own_ranks.exists()), in which case
   ownership is decided to favor ranks which
   have fewer objects, with smallest rank being
   the tiebreaker when ranks have the same number
   of objects.
   if (own_ranks) is specified, it dictates the owner
   ranks for all copies and is assumed to be consistent
   over duplicate copies.
   in this case we use that information to compute the
   owner indices for all copies.

   in both cases, the ownership is computed
   by linearly partitioning the global numbers
   across the communicator processes, transmitting
   information to the linearly partitioned global
   objects, and deciding ownership in this partitioning.
   this can be viewed as bootstrapping the true ownership
   from a more crude temporary ownership based on
   linear partitioning. */

Remotes owners_from_globals(CommPtr comm,
    Read<GO> globals, Read<I32> own_ranks);
