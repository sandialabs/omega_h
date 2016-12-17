#ifndef ARRAY_HPP
#define ARRAY_HPP

#include "Omega_h_functors.hpp"
#include "internal.hpp"

namespace Omega_h {

template <class T>
bool operator==(Read<T> a, Read<T> b);

template <class T>
Write<T> deep_copy(Read<T> a);

template <typename T>
typename StandinTraits<T>::type get_sum(Read<T> a);
template <typename T>
T get_min(Read<T> a);
template <typename T>
T get_max(Read<T> a);

template <typename T>
typename StandinTraits<T>::type get_sum(CommPtr comm, Read<T> a);
template <typename T>
T get_min(CommPtr comm, Read<T> a);
template <typename T>
T get_max(CommPtr comm, Read<T> a);

template <typename T>
struct MinMax {
  T min;
  T max;
};

template <typename T>
MinMax<T> get_minmax(CommPtr comm, Read<T> a);

bool are_close(Reals a, Reals b, Real tol = EPSILON, Real floor = EPSILON);

/* "a" may be larger than "b" by some integer factor */
template <typename T>
Read<T> multiply_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> divide_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> subtract_each(Read<T> a, Read<T> b);
template <typename T>
Read<T> add_to_each(Read<T> a, T b);
template <typename T>
Read<I8> each_geq_to(Read<T> a, T b);
template <typename T>
Read<I8> each_gt(Read<T> a, T b);
template <typename T>
Read<I8> each_lt(Read<T> a, T b);
template <typename T>
Read<I8> gt_each(Read<T> a, Read<T> b);
template <typename T>
Read<I8> each_neq_to(Read<T> a, T b);
Read<I8> land_each(Read<I8> a, Read<I8> b);
Read<I8> lor_each(Read<I8> a, Read<I8> b);
Read<I8> bit_or_each(Read<I8> a, Read<I8> b);
Read<I8> bit_neg_each(Read<I8> a);

template <typename T>
Read<T> get_component(Read<T> a, Int ncomps, Int comp);

#define INST_DECL(T)                                                           \
  extern template class Write<T>;                                              \
  extern template class Read<T>;                                               \
  extern template class HostWrite<T>;                                          \
  extern template class HostRead<T>;                                           \
  extern template bool operator==(Read<T> a, Read<T> b);                       \
  extern template typename StandinTraits<T>::type get_sum(Read<T> a);              \
  extern template T get_min(Read<T> a);                                            \
  extern template T get_max(Read<T> a);                                            \
  extern template typename StandinTraits<T>::type get_sum(                         \
      CommPtr comm, Read<T> a);                                                \
  extern template T get_min(CommPtr comm, Read<T> a);                              \
  extern template T get_max(CommPtr comm, Read<T> a);                              \
  extern template MinMax<T> get_minmax(CommPtr comm, Read<T> a);                              \
  extern template Write<T> deep_copy(Read<T> a);                               \
  extern template Read<T> multiply_each(Read<T> a, Read<T> b);                 \
  extern template Read<T> divide_each(Read<T> a, Read<T> b);                   \
  extern template Read<T> add_each(Read<T> a, Read<T> b);                      \
  extern template Read<T> subtract_each(Read<T> a, Read<T> b);                 \
  extern template Read<T> add_to_each(Read<T> a, T b);                         \
  extern template Read<I8> each_geq_to(Read<T> a, T b);                        \
  extern template Read<I8> each_gt(Read<T> a, T b);                            \
  extern template Read<I8> each_lt(Read<T> a, T b);                            \
  extern template Read<I8> each_neq_to(Read<T> a, T b);                        \
  extern template Read<I8> gt_each(Read<T> a, Read<T> b);                      \
  extern template Read<T> get_component(Read<T> a, Int ncomps, Int comp);

INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL

}  // end namespace Omega_h

#endif
