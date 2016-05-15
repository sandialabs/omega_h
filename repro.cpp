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
  double unit_;
  ReproSum(Reals a, double unit):a_(a),unit_(unit) {}
  INLINE void operator()(UInt i, value_type& update) const {
    update = update + Int128::from_double(a_[i], unit_);
  }
};

Real repro_sum(Reals a, int expo) {
  double unit = exp2(double(expo - MANTISSA_BITS));
  return parallel_reduce(a.size(), ReproSum(a, unit)).to_double(unit);
}

Real repro_sum(Reals a) {
  return repro_sum(a, max_exponent(a));
}
