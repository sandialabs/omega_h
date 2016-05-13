struct MaxExponent : public MaxFunctor<int> {
  Reals a_;
  MaxExponent(Reals a):a_(a) {}
  INLINE void operator()(UInt i, value_type& update) const {
    int expo;
    frexp(a_[i], &expo);
    if (expo > update)
      update = expo;
  }
};

int max_exponent(Reals a) {
  return parallel_reduce(a.size(), MaxExponent(a));
}

struct ReproSum : public SumFunctor<Int128> {
  Reals a_;
  int expo_;
  ReproSum(Reals a, int expo):a_(a),expo_(expo) {}
  INLINE void operator()(UInt i, value_type& update) const {
    update = update + Int128(a_[i], expo_);
  }
};

Real repro_sum(Reals a) {
  int expo = max_exponent(a);
  return parallel_reduce(a.size(), ReproSum(a, expo)).as_double(expo);
}
