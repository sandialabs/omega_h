#ifndef OMEGA_H_FUNCTORS_HPP
#define OMEGA_H_FUNCTORS_HPP

#include "Omega_h_math.hpp"

namespace Omega_h {

/* values smaller than 64 bits cause wrong reduction
 * answers on GPUs; the StandinTraits system raises
 * such values to their equivalent 64-bit type which is
 * then used as the reduction value_type
 * in addition, we use this type for global MPI reductions,
 * so it makes even more sense for it to be 64 bits.
 */
template <typename T>
struct StandinTraits {
  typedef T type;
};

template <>
struct StandinTraits<I8> {
  typedef I64 type;
};

template <>
struct StandinTraits<I32> {
  typedef I64 type;
};

struct AndFunctor {
  typedef I64 value_type;
  OMEGA_H_INLINE void init(value_type& update) const { update = 1; }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = update && input;
  }
};

template <typename T>
struct MaxFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  OMEGA_H_INLINE void init(value_type& update) const {
    update = ArithTraits<T>::min();
  }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = max2(update, input);
  }
};

template <typename T>
struct MinFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  OMEGA_H_INLINE void init(value_type& update) const {
    update = ArithTraits<T>::max();
  }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = min2(update, input);
  }
};

template <typename T>
struct SumFunctor {
  using value_type = typename StandinTraits<T>::type;
  using input_type = T;
  OMEGA_H_INLINE void init(value_type& update) const { update = 0; }
  OMEGA_H_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = update + input;
  }
};

}  // end namespace Omega_h

#endif
