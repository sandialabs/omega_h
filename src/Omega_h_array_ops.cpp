#include "Omega_h_array_ops.hpp"

#include "Omega_h_internal.hpp"

namespace Omega_h {

template <class T>
struct SameContent : public AndFunctor {
  Read<T> a_;
  Read<T> b_;
  SameContent(Read<T> a, Read<T> b) : a_(a), b_(b) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = update && (a_[i] == b_[i]);
  }
};

template <class T>
bool operator==(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  return parallel_reduce(a.size(), SameContent<T>(a, b));
}

template <typename T>
struct Sum : public SumFunctor<T> {
  using typename SumFunctor<T>::value_type;
  Read<T> a_;
  Sum(Read<T> a) : a_(a) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = update + a_[i];
  }
};

template <typename T>
typename StandinTraits<T>::type get_sum(Read<T> a) {
  return parallel_reduce(a.size(), Sum<T>(a));
}

template <typename T>
typename StandinTraits<T>::type get_sum(CommPtr comm, Read<T> a) {
  return comm->allreduce(get_sum(a), OMEGA_H_SUM);
}

template <typename T>
struct Min : public MinFunctor<T> {
  using typename MinFunctor<T>::value_type;
  Read<T> a_;
  Min(Read<T> a) : a_(a) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = min2<value_type>(update, a_[i]);
  }
};

template <typename T>
T get_min(Read<T> a) {
  auto r = parallel_reduce(a.size(), Min<T>(a));
  return static_cast<T>(r);  // see StandinTraits
}

template <typename T>
struct Max : public MaxFunctor<T> {
  using typename MaxFunctor<T>::value_type;
  Read<T> a_;
  Max(Read<T> a) : a_(a) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = max2<value_type>(update, a_[i]);
  }
};

template <typename T>
T get_max(Read<T> a) {
  auto r = parallel_reduce(a.size(), Max<T>(a));
  return static_cast<T>(r);  // see StandinTraits
}

template <typename T>
T get_min(CommPtr comm, Read<T> a) {
  return comm->allreduce(get_min(a), OMEGA_H_MIN);
}

template <typename T>
T get_max(CommPtr comm, Read<T> a) {
  return comm->allreduce(get_max(a), OMEGA_H_MAX);
}

template <typename T>
MinMax<T> get_minmax(CommPtr comm, Read<T> a) {
  return {get_min(comm, a), get_max(comm, a)};
}

struct AreClose : public AndFunctor {
  Reals a_;
  Reals b_;
  Real tol_;
  Real floor_;
  AreClose(Reals a, Reals b, Real tol, Real floor)
      : a_(a), b_(b), tol_(tol), floor_(floor) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = update && are_close(a_[i], b_[i], tol_, floor_);
  }
};

bool are_close(Reals a, Reals b, Real tol, Real floor) {
  CHECK(a.size() == b.size());
  return static_cast<bool>(
      parallel_reduce(a.size(), AreClose(a, b, tol, floor)));
}

struct AreCloseAbs : public AndFunctor {
  Reals a_;
  Reals b_;
  Real tol_;
  AreCloseAbs(Reals a, Reals b, Real tol) : a_(a), b_(b), tol_(tol) {}
  DEVICE void operator()(LO i, value_type& update) const {
    update = update && (fabs(a_[i] - b_[i]) <= tol_);
  }
};

bool are_close_abs(Reals a, Reals b, Real tol) {
  CHECK(a.size() == b.size());
  return static_cast<bool>(parallel_reduce(a.size(), AreCloseAbs(a, b, tol)));
}

template <typename T>
Read<T> multiply_each_by(T factor, Read<T> a) {
  Write<T> b(a.size());
  auto f = LAMBDA(LO i) { b[i] = a[i] * factor; };
  parallel_for(a.size(), f);
  return b;
}

template <typename T>
Read<T> multiply_each(Read<T> a, Read<T> b) {
  if (b.size() == 0) {
    CHECK(a.size() == 0);
    return a;
  }
  auto width = divide_no_remainder(a.size(), b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) {
    for (Int j = 0; j < width; ++j) {
      c[i * width + j] = a[i * width + j] * b[i];
    }
  };
  parallel_for(b.size(), f);
  return c;
}

template <typename T>
Read<T> divide_each(Read<T> a, Read<T> b) {
  if (b.size() == 0) {
    CHECK(a.size() == 0);
    return a;
  }
  auto width = divide_no_remainder(a.size(), b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) {
    for (Int j = 0; j < width; ++j) {
      c[i * width + j] = a[i * width + j] / b[i];
    }
  };
  parallel_for(b.size(), f);
  return c;
}

template <typename T>
Read<T> divide_each_by(T factor, Read<T> a) {
  Write<T> b(a.size());
  auto f = LAMBDA(LO i) { b[i] = a[i] / factor; };
  parallel_for(a.size(), f);
  return b;
}

template <typename T>
Read<T> add_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = a[i] + b[i]; };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> subtract_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = a[i] - b[i]; };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> add_to_each(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = a[i] + b; };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> subtract_from_each(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = a[i] - b; };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_geq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] >= b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_leq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] <= b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_gt(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] > b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_lt(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] < b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> gt_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] > b[i]); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> geq_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] >= b[i]); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> min_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = min2(a[i], b[i]); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> max_each(Read<T> a, Read<T> b) {
  CHECK(a.size() == b.size());
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = max2(a[i], b[i]); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<T> each_max_with(Read<T> a, T b) {
  Write<T> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = max2(a[i], b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_neq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] != b); };
  parallel_for(c.size(), f);
  return c;
}

template <typename T>
Read<I8> each_eq_to(Read<T> a, T b) {
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] == b); };
  parallel_for(c.size(), f);
  return c;
}

Read<I8> land_each(Read<I8> a, Read<I8> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] && b[i]); };
  parallel_for(c.size(), f);
  return c;
}

Read<I8> lor_each(Read<I8> a, Read<I8> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] || b[i]); };
  parallel_for(c.size(), f);
  return c;
}

Read<I8> bit_or_each(Read<I8> a, Read<I8> b) {
  CHECK(a.size() == b.size());
  Write<I8> c(a.size());
  auto f = LAMBDA(LO i) { c[i] = (a[i] | b[i]); };
  parallel_for(c.size(), f);
  return c;
}

Read<I8> bit_neg_each(Read<I8> a) {
  Write<I8> b(a.size());
  auto f = LAMBDA(LO i) { b[i] = ~(a[i]); };
  parallel_for(a.size(), f);
  return b;
}

Read<Real> fabs_each(Read<Real> a) {
  Write<I8> b(a.size());
  auto f = LAMBDA(LO i) { b[i] = fabs(a[i]); };
  parallel_for(a.size(), f);
  return b;
}

template <typename T>
Read<T> get_component(Read<T> a, Int ncomps, Int comp) {
  Write<T> b(divide_no_remainder(a.size(), ncomps));
  auto f = LAMBDA(LO i) { b[i] = a[i * ncomps + comp]; };
  parallel_for(b.size(), f);
  return b;
}

template <typename T>
void set_component(Write<T> out, Read<T> a, Int ncomps, Int comp) {
  auto f = LAMBDA(LO i) { out[i * ncomps + comp] = a[i]; };
  parallel_for(a.size(), f);
}

template <typename T>
struct FindLast : public MaxFunctor<LO> {
  using typename MaxFunctor<LO>::value_type;
  Read<T> array_;
  T value_;
  FindLast(Read<T> array, T value) : array_(array), value_(value) {}
  DEVICE void operator()(LO i, value_type& update) const {
    if (array_[i] == value_) {
      update = max2<value_type>(update, i);
    }
  }
};

template <typename T>
LO find_last(Read<T> array, T value) {
  return static_cast<LO>(
      parallel_reduce(array.size(), FindLast<T>(array, value)));
}

/* A reproducible sum of floating-point values.
   this operation is one of the key places where
   a program's output begins to depend on parallel
   partitioning and traversal order, because
   floating-point values do not produce the same
   sum when added in a different order.

   IEEE 754 64-bit floating point format is assumed,
   which has 52 bits in the fraction.

   The idea here is to add the numbers as fixed-point values.
   max_exponent() finds the largest exponent (e) such that
   all values are (<= 2^(e)).
   We then use the value (2^(e - 52)) as the unit, and sum all
   values as integers in that unit.
   This is guaranteed to be at least as accurate as the
   worst-case ordering of the values, i.e. being added
   in order of decreasing magnitude.

   If we used a 64-bit integer type, we would only be
   able to reliably add up to (2^12 = 4096) values
   (64 - 52 = 12).
   Thus we use a 128-bit integer type.
   This allows us to reliably add up to (2^76 > 10^22)
   values. By comparison, supercomputers today
   support a maximum of one million MPI ranks (10^6)
   and each rank typically can't hold more than
   one billion values (10^9), for a total of (10^15) values.
*/

struct MaxExponent : public MaxFunctor<int> {
  Reals a_;
  MaxExponent(Reals a) : a_(a) {}
  DEVICE void operator()(Int i, value_type& update) const {
    int expo;
    frexp(a_[i], &expo);
    if (expo > update) update = expo;
  }
};

static int max_exponent(Reals a) {
  return static_cast<int>(parallel_reduce(a.size(), MaxExponent(a)));
}

struct ReproSum : public SumFunctor<Int128> {
  Reals a_;
  double unit_;
  ReproSum(Reals a, double unit) : a_(a), unit_(unit) {}
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
  int expo = comm->allreduce(max_exponent(a), OMEGA_H_MAX);
  double unit = exp2(double(expo - MANTISSA_BITS));
  Int128 fixpt_sum = parallel_reduce(a.size(), ReproSum(a, unit));
  fixpt_sum = comm->add_int128(fixpt_sum);
  return fixpt_sum.to_double(unit);
}

void repro_sum(CommPtr comm, Reals a, Int ncomps, Real result[]) {
  for (Int comp = 0; comp < ncomps; ++comp) {
    result[comp] = repro_sum(comm, get_component(a, ncomps, comp));
  }
}

Reals interpolate_between(Reals a, Reals b, Real t) {
  CHECK(a.size() == b.size());
  auto n = a.size();
  auto out = Write<Real>(n);
  auto f = LAMBDA(LO i) { out[i] = a[i] * (1.0 - t) + b[i] * t; };
  parallel_for(n, f);
  return out;
}

template <typename Tout, typename Tin>
Read<Tout> array_cast(Read<Tin> in) {
  auto out = Write<Tout>(in.size());
  auto f = LAMBDA(LO i) { out[i] = static_cast<Tout>(in[i]); };
  parallel_for(in.size(), f);
  return out;
}

#define INST(T)                                                                \
  template bool operator==(Read<T> a, Read<T> b);                              \
  template typename StandinTraits<T>::type get_sum(Read<T> a);                 \
  template T get_min(Read<T> a);                                               \
  template T get_max(Read<T> a);                                               \
  template typename StandinTraits<T>::type get_sum(CommPtr comm, Read<T> a);   \
  template T get_min(CommPtr comm, Read<T> a);                                 \
  template T get_max(CommPtr comm, Read<T> a);                                 \
  template MinMax<T> get_minmax(CommPtr comm, Read<T> a);                      \
  template Read<T> multiply_each_by(T factor, Read<T> x);                      \
  template Read<T> divide_each_by(T factor, Read<T> x);                        \
  template Read<T> multiply_each(Read<T> a, Read<T> b);                        \
  template Read<T> divide_each(Read<T> a, Read<T> b);                          \
  template Read<T> add_each(Read<T> a, Read<T> b);                             \
  template Read<T> subtract_each(Read<T> a, Read<T> b);                        \
  template Read<T> min_each(Read<T> a, Read<T> b);                             \
  template Read<T> max_each(Read<T> a, Read<T> b);                             \
  template Read<T> each_max_with(Read<T> a, T b);                              \
  template Read<T> add_to_each(Read<T> a, T b);                                \
  template Read<T> subtract_from_each(Read<T> a, T b);                         \
  template Read<I8> each_geq_to(Read<T> a, T b);                               \
  template Read<I8> each_leq_to(Read<T> a, T b);                               \
  template Read<I8> each_gt(Read<T> a, T b);                                   \
  template Read<I8> each_lt(Read<T> a, T b);                                   \
  template Read<I8> each_neq_to(Read<T> a, T b);                               \
  template Read<I8> each_eq_to(Read<T> a, T b);                                \
  template Read<I8> gt_each(Read<T> a, Read<T> b);                             \
  template Read<T> get_component(Read<T> a, Int ncomps, Int comp);             \
  template void set_component(Write<T> out, Read<T> a, Int ncomps, Int comp);             \
  template LO find_last(Read<T> array, T value);

INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

template Read<Real> array_cast(Read<I32>);
template Read<I32> array_cast(Read<I8>);

}  // end namespace Omega_h
