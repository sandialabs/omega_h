struct MaxExponent : public MaxFunctor<int> {
  Reals a_;
  MaxExponent(Reals a):a_(a) {}
  DEVICE void operator()(Int i, value_type& update) const {
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
  DEVICE void operator()(Int i, value_type& update) const {
    update = update + Int128::from_double(a_[i], unit_);
  }
};

Real repro_sum(Reals a) {
  int expo = max_exponent(a);
  double unit = exp2(double(expo - MANTISSA_BITS));
  Int128 fixpt_sum = parallel_reduce(a.size(), ReproSum(a, unit));
  return fixpt_sum.to_double(unit);
}

Real repro_sum(CommPtr comm, Reals a) {
  int expo = comm->allreduce(max_exponent(a), MAX);
  double unit = exp2(double(expo - MANTISSA_BITS));
  Int128 fixpt_sum = parallel_reduce(a.size(), ReproSum(a, unit));
  fixpt_sum = comm->add_int128(fixpt_sum);
  return fixpt_sum.to_double(unit);
}

void repro_sum(CommPtr comm, Reals a, Int ncomps, Real result[]) {
  CHECK(a.size() % ncomps == 0);
  auto n = a.size() / ncomps;
  for (Int j = 0; j < ncomps; ++j) {
    Write<Real> comp(n);
    auto f = LAMBDA(LO i) {
      comp[i] = a[i * ncomps + j];
    };
    parallel_for(n, f);
    result[j] = repro_sum(comm, Reals(comp));
  }
}
