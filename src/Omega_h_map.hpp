#ifndef OMEGA_H_MAP_HPP
#define OMEGA_H_MAP_HPP

#include <Omega_h_array.hpp>
#include <Omega_h_graph.hpp>

namespace Omega_h {

template <typename T>
void add_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width);

template <typename T>
void map_into(Read<T> a_data, LOs a2b, Write<T> b_data, Int width);

template <typename T>
Read<T> map_onto(Read<T> a_data, LOs a2b, LO nb, T init_val, Int width);

template <typename T>
Read<T> unmap(LOs a2b, Read<T> b_data, Int width);

template <typename T>
Read<T> expand(Read<T> a_data, LOs a2b, Int width);

template <typename T>
Read<T> permute(Read<T> a_data, LOs a2b, Int width);

LOs multiply_fans(LOs a2b, LOs a2c);

LOs compound_maps(LOs a2b, LOs b2c);

LOs invert_permutation(LOs a2b);

Read<I8> invert_marks(Read<I8> marks);

LOs collect_marked(Read<I8> marks);

Read<I8> mark_image(LOs a2b, LO nb);

/* The map (a2b) is injective iff each value of a maps to (a) unique value of (b).
   Let b2a = invert_injective_map(a2b, nb);
   Then b2a[a2b[a]] == a for all a in [0, a2b.size() - 1]
   For values of (b) that don't have an inverse, b2a[b] == -1
 */
LOs invert_injective_map(LOs a2b, LO nb);

LOs invert_funnel(LOs ab2a, LO na);

/* Inverts a map which is not necessarily injective.
   Since multiple (a)s may map to the same (b), the
   result is a "graph" structure.
   Note that the map (a2b) is essentially a graph where
   all nodes in the A set have exactly one outgoing edge.
   The inversion is done by sorting these "edges"
   (the entries of (a2b)) by their destination node (b).
   Because we use a stable sort, in the resulting graph
   the edges leaving a single (b) are sorted by their
   source node (a), because that was their relative
   order in (a2b) */
Graph invert_map(LOs a2b, LO nb);

LOs get_degrees(LOs offsets);

LOs invert_fan(LOs a2b);

template <typename T>
Read<T> fan_sum(LOs a2b, Read<T> b_data);
template <typename T>
Read<T> fan_max(LOs a2b, Read<T> b_data);
template <typename T>
Read<T> fan_min(LOs a2b, Read<T> b_data);
template <typename T>
Read<T> fan_reduce(LOs a2b, Read<T> b_data, Int width, Omega_h_Op op);

#define INST_T(T)                                                              \
  extern template Read<T> permute(Read<T> a_data, LOs a2b, Int width);         \
  extern template void add_into(                                               \
      Read<T> a_data, LOs a2b, Write<T> b_data, Int width);                    \
  extern template void map_into(                                               \
      Read<T> a_data, LOs a2b, Write<T> b_data, Int width);                    \
  extern template Read<T> map_onto(                                            \
      Read<T> a_data, LOs a2b, LO nb, T, Int width);                           \
  extern template Read<T> unmap(LOs a2b, Read<T> b_data, Int width);           \
  extern template Read<T> expand(Read<T> a_data, LOs a2b, Int width);          \
  extern template Read<T> fan_reduce(                                          \
      LOs a2b, Read<T> b_data, Int width, Omega_h_Op op);
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T

}  // end namespace Omega_h

#endif
