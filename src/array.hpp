#ifndef ARRAY_HPP
#define ARRAY_HPP

#include "omega_h_functors.hpp"
#include "internal.hpp"

namespace osh {

template <class T>
bool operator==(Read<T> a, Read<T> b);

template <class T>
Write<T> deep_copy(Read<T> a);

template <typename T>
typename StandinTraits<T>::type sum(Read<T> a);
template <typename T>
T min(Read<T> a);
template <typename T>
T max(Read<T> a);

template <typename T>
typename StandinTraits<T>::type sum(CommPtr comm, Read<T> a);
template <typename T>
T min(CommPtr comm, Read<T> a);
template <typename T>
T max(CommPtr comm, Read<T> a);

bool are_close(Reals a, Reals b, Real tol = EPSILON, Real floor = EPSILON);

template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a);
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
template <typename T>
Read<I8> each_eq_to(Read<T> a, T b);
Read<I8> land_each(Read<I8> a, Read<I8> b);
Read<I8> lor_each(Read<I8> a, Read<I8> b);

template <typename T>
Read<T> get_component(Read<T> a, Int ncomps, Int comp);

#define INST_DECL(T)                                                \
  extern template class Write<T>;                                   \
  extern template class Read<T>;                                    \
  extern template class HostWrite<T>;                               \
  extern template class HostRead<T>;                                \
  extern template bool operator==(Read<T> a, Read<T> b);            \
  extern template typename StandinTraits<T>::type sum(Read<T> a);   \
  extern template T min(Read<T> a);                                 \
  extern template T max(Read<T> a);                                 \
  extern template typename StandinTraits<T>::type sum(CommPtr comm, \
                                                      Read<T> a);   \
  extern template T min(CommPtr comm, Read<T> a);                   \
  extern template T max(CommPtr comm, Read<T> a);                   \
  extern template Write<T> deep_copy(Read<T> a);                    \
  extern template Read<T> multiply_each_by(T factor, Read<T> x);    \
  extern template Read<T> multiply_each(Read<T> a, Read<T> b);      \
  extern template Read<T> divide_each(Read<T> a, Read<T> b);        \
  extern template Read<T> add_each(Read<T> a, Read<T> b);           \
  extern template Read<T> subtract_each(Read<T> a, Read<T> b);      \
  extern template Read<T> add_to_each(Read<T> a, T b);              \
  extern template Read<I8> each_geq_to(Read<T> a, T b);             \
  extern template Read<I8> each_gt(Read<T> a, T b);                 \
  extern template Read<I8> each_lt(Read<T> a, T b);                 \
  extern template Read<I8> each_neq_to(Read<T> a, T b);             \
  extern template Read<I8> each_eq_to(Read<T> a, T b);              \
  extern template Read<I8> gt_each(Read<T> a, Read<T> b);           \
  extern template Read<T> get_component(Read<T> a, Int ncomps, Int comp);

INST_DECL(I8)
INST_DECL(I32)
INST_DECL(I64)
INST_DECL(Real)
#undef INST_DECL

}  // end namespace osh

#endif
