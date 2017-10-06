#ifndef OMEGA_H_SORT_HPP
#define OMEGA_H_SORT_HPP

#include <Omega_h_array.hpp>

namespace Omega_h {

/* Compute the permutation which sorts the given keys.
   Each key is a tuple of N integers of type T.
   T may be 32 or 64 bits, and N may be 1, 2, or 3.
   Let perm = sort_by_keys(keys, width);
   Then the key at perm[i] is <= the key at perm[i + 1].
   In other words, old_key_index = perm[new_key_index]
   The sorting algorithm used is stable, so if two keys
   are equal then their relative order after permutation
   will be the same as before the permutation.
   Tuples (keys) are sorted into lexical order, so they
   will be sorted by the first integer first, second
   integer second, etc.
 */
template <typename T>
LOs sort_by_keys(Read<T> keys, Int width = 1);

#define OMEGA_H_INST_DECL(T)                                                   \
  extern template LOs sort_by_keys(Read<T> keys, Int width);
OMEGA_H_INST_DECL(LO)
OMEGA_H_INST_DECL(GO)
#undef OMEGA_H_INST_DECL

template <typename T>
void sort_small_range(
    Read<T> items2values, LOs* p_perm, LOs* p_fan, Read<T>* p_uniq);

extern template void sort_small_range(
    Read<I32> items2values, LOs* p_perm, LOs* p_fan, Read<I32>* p_uniq);

}  // end namespace Omega_h

#endif
