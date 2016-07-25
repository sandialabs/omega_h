#ifndef FUNCTORS_HPP
#define FUNCTORS_HPP

#include "algorithm.hpp"
#include "traits.hpp"

namespace osh {

/* we use I8 for storing very small values efficiently,
   but since it boils down to 'char' usually, this causes
   problems when printing values or using Kokkos reductions.
   StandinTraits just replaces I8 with I32 */

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
  INLINE void init(value_type& update) const { update = 1; }
  INLINE void join(volatile value_type& update,
                   const volatile value_type& input) const {
    update = update && input;
  }
};

template <typename T>
struct MaxFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  INLINE void init(value_type& update) const { update = ArithTraits<T>::min(); }
  INLINE void join(volatile value_type& update,
                   const volatile value_type& input) const {
    update = max2(update, input);
  }
};

template <typename T>
struct MinFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  INLINE void init(value_type& update) const { update = ArithTraits<T>::max(); }
  INLINE void join(volatile value_type& update,
                   const volatile value_type& input) const {
    update = min2(update, input);
  }
};

template <typename T>
struct SumFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  INLINE void init(value_type& update) const { update = 0; }
  INLINE void join(volatile value_type& update,
                   const volatile value_type& input) const {
    update = update + input;
  }
};

}  // end namespace osh

#endif
