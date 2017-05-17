#ifndef OMEGA_H_BIPART_HPP
#define OMEGA_H_BIPART_HPP

#include <Omega_h_dist.hpp>

namespace Omega_h {

/* given a distributed set of items, and a boolean
   for each item indicating which half of the new
   partitioning it belongs in, construct a Dist
   object that maps current items to their destinations,
   such that the first half of the ranks have
   all the items marked 0 and the second half
   of the ranks have all the items marked 1.
   this is a useful subroutine for parallel
   sort-like operations (RIB in particular) */

Dist bi_partition(CommPtr comm, Read<I8> marks);

}  // end namespace Omega_h

#endif
