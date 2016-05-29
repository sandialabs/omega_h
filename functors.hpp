struct AndFunctor {
  typedef Int value_type;
  INLINE void init(value_type& update) const {
    update = 1;
  }
  INLINE void join(volatile value_type& update,
      const volatile value_type& input) const {
    update = update && input;
  }
};

template <typename T>
struct MaxFunctor {
  typedef T value_type;
  INLINE void init(T& update) const {
    update = ArithTraits<T>::min();
  }
  INLINE void join(volatile T& update, const volatile T& input) const {
    update = max2(update, input);
  }
};

template <typename T>
struct MinFunctor {
  typedef T value_type;
  INLINE void init(T& update) const {
    update = ArithTraits<T>::max();
  }
  INLINE void join(volatile T& update, const volatile T& input) const {
    update = min2(update, input);
  }
};

template <typename T>
struct SumFunctor {
  typedef T value_type;
  INLINE void init(T& update) const {
    update = 0;
  }
  INLINE void join(volatile T& update, const volatile T& input) const {
    update = update + input;
  }
};
