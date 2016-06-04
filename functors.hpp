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
  typedef I32 type;
};

struct AndFunctor {
  typedef Int value_type;
  OSH_INLINE void init(value_type& update) const {
    update = 1;
  }
  OSH_INLINE void join(volatile value_type& update,
      const volatile value_type& input) const {
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
      volatile value_type& update,
      const volatile value_type& input) const {
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
      volatile value_type& update,
      const volatile value_type& input) const {
    update = min2(update, input);
  }
};

template <typename T>
struct SumFunctor {
  typedef typename StandinTraits<T>::type value_type;
  typedef T input_type;
  OSH_INLINE void init(value_type& update) const {
    update = 0;
  }
  OSH_INLINE void join(
      volatile value_type& update,
      const volatile value_type& input) const {
    update = update + input;
  }
};
