#ifndef OMEGA_H_FUNCTORS_HPP
#define OMEGA_H_FUNCTORS_HPP

#include "omega_h_math.hpp"

namespace osh {

/* values smaller than 64 bits cause wrong reduction
 * answers on GPUs; the StandinTraits system raises
 * such values to their equivalent 64-bit type which is
 * then used as the reduction value_type
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
  OSH_INLINE void init(value_type& update) const { update = 1; }
  OSH_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = update && input;
  }
};

template <typename T>
struct MaxFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  OSH_INLINE void init(value_type& update) const {
    update = ArithTraits<T>::min();
  }
  OSH_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = max2(update, input);
  }
};

template <typename T>
struct MinFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  OSH_INLINE void init(value_type& update) const {
    update = ArithTraits<T>::max();
  }
  OSH_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = min2(update, input);
  }
};

template <typename T>
struct SumFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  OSH_INLINE void init(value_type& update) const { update = 0; }
  OSH_INLINE void join(
      volatile value_type& update, const volatile value_type& input) const {
    update = update + input;
  }
};

}  // end namespace osh

#endif
