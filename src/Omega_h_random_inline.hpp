#ifndef OMEGA_H_RANDOM_INLINE_HPP
#define OMEGA_H_RANDOM_INLINE_HPP

#include <random123/philox.h>
#include <Omega_h_few.hpp>
#include <Omega_h_scalar.hpp>
#include <Omega_h_vector.hpp>

namespace Omega_h {

OMEGA_H_INLINE Few<std::uint64_t, 2> run_philox_cbrng(
    Few<std::uint64_t, 2> ctr, std::uint64_t key) {
  using PhiloxType = r123::Philox4x32;
  using ctr_type = typename PhiloxType::ctr_type;
  using key_type = typename PhiloxType::key_type;
  PhiloxType philox_cbrng;
  ctr_type ctr_philox;
  key_type key_philox;
  ctr_philox[0] = static_cast<std::uint32_t>(ctr[0]);
  ctr_philox[1] = static_cast<std::uint32_t>(ctr[0] >> 32);
  ctr_philox[2] = static_cast<std::uint32_t>(ctr[1]);
  ctr_philox[3] = static_cast<std::uint32_t>(ctr[1] >> 32);
  key_philox[0] = static_cast<std::uint32_t>(key);
  key_philox[1] = static_cast<std::uint32_t>(key >> 32);
  ctr_philox = philox_cbrng(ctr_philox, key_philox);
  ctr[0] = ctr_philox[1];
  ctr[0] <<= 32;
  ctr[0] |= ctr_philox[0];
  ctr[1] = ctr_philox[3];
  ctr[1] <<= 32;
  ctr[1] |= ctr_philox[2];
  return ctr;
}

OMEGA_H_INLINE Real unit_uniform_deviate_from_uint64(std::uint64_t x) {
  return static_cast<Real>(x) /
         static_cast<Real>(ArithTraits<std::uint64_t>::max());
}

class UnitUniformDistribution {
 public:
  OMEGA_H_INLINE UnitUniformDistribution(
      I64 seed_in, I64 key_in, I64 counter_in) {
    counter[1] = static_cast<std::uint64_t>(seed_in);
    key = static_cast<std::uint64_t>(key_in);
    counter[0] = static_cast<std::uint64_t>(counter_in);
    if (counter[0] % 2)
      even_step();
    else
      state[1] = 0;  // Silences GCC 7.2.0 -Wmaybe-uninitialized about state[1]
                     // in operator()
  }
  OMEGA_H_INLINE Real operator()() {
    Real ret;
    if (counter[0] % 2) {
      ret = unit_uniform_deviate_from_uint64(state[1]);
    } else {
      even_step();
      ret = unit_uniform_deviate_from_uint64(state[0]);
    }
    ++counter[0];
    return ret;
  }

 private:
  OMEGA_H_INLINE void even_step() {
    counter[0] >>= 1;
    state = run_philox_cbrng(counter, key);
    counter[0] <<= 1;
  }
  Few<std::uint64_t, 2> counter;
  std::uint64_t key;
  Few<std::uint64_t, 2> state;
};

/* I was initially inspired by GCC's implementation of std::normal_distribution,
   which is in turn taken from:

   Devroye, Luc. "Non-Uniform Random Variate Generation".
   Springer, New York, NY, 1986. 1-26.
   Chapter 4: "Non-Uniform Random Variate Generation"
   Section V.4: "Polar Method"
   SubSection 4.4: "Generating normal random variates in batches"
   Page 235

   However, the rejection method used to sample the unit circle would mean
   that counters advance differently for different random streams, which would
   complicate programming.

   Thus, we return to the older but more closed-form approach from:
   G. E. P. Box and Mervin E. Muller
   "A Note on the Generation of Random Normal Deviates"
   The Annals of Mathematical Statistics, 1958

 */

class StandardNormalDistribution {
 public:
  OMEGA_H_INLINE StandardNormalDistribution() : counter(0) {}
  OMEGA_H_INLINE Real operator()(UnitUniformDistribution& uniform_rng) {
    Real ret;
    if (counter) {
      ret = state[1];
    } else {
      auto U_1 = uniform_rng();
      auto U_2 = uniform_rng();
      /* this is the suspicious part. we can't feed exact zero to std::log,
         or we'll get a pole error. other implementations will use rejection
         here, but again I want the number of uniform deviates generated
         to be a known constant, not a function of their values */
      if (U_1 == 0.0) U_1 = DBL_MIN;
      auto common = std::sqrt(-2.0 * std::log(U_1));
      state[0] = common * std::cos(2.0 * Omega_h::PI * U_2);
      state[1] = common * std::sin(2.0 * Omega_h::PI * U_2);
      ret = state[0];
    }
    return ret;
  }

 private:
  int counter;
  Vector<2> state;
};

/* fortunately, the Weibull distribution has a closed-form inverse CDF
(quantile)

https://en.wikipedia.org/wiki/Weibull_distribution#Cumulative_distribution_function
 */
OMEGA_H_INLINE Real weibull_quantile(Real shape, Real scale, Real p) {
  auto one_minus_p = 1.0 - p;
  if (one_minus_p == 0) one_minus_p = DBL_MIN;
  return scale * std::pow(double(-std::log(one_minus_p)), double(1.0 / shape));
}

OMEGA_H_INLINE Real standard_normal_density(Real value) {
  constexpr auto one_over_sqrt_two_pi = 0.3989422804014327;
  return one_over_sqrt_two_pi * std::exp((-1.0 / 2.0) * square(value));
}

OMEGA_H_INLINE Real general_normal_density(
    Real mean, Real standard_deviation, Real value) {
  return (1.0 / standard_deviation) *
         standard_normal_density((value - mean) / standard_deviation);
}

OMEGA_H_INLINE Real weibull_density(Real shape, Real scale, Real value) {
  if (value < 0.0) return 0.0;
  return (shape / scale) * std::pow(value / scale, shape - 1.0) *
         std::exp(-std::pow(value / scale, shape));
}

// regularized lower incomplete gamma function, by series expansion
OMEGA_H_INLINE Real regularized_lower_incomplete_gamma(Real s, Real z) {
  Real sum = 1.0;
  Real x = 1.0;
  for (Int k = 1; k < 100; ++k) {
    x *= z / (s + k);
    sum += x;
    if (sum > 1e100) return 1.0;
    if (x / sum < 1e-14) break;
  }
  auto a = s * std::log(z) - z - std::lgamma(s + 1.0);
  return std::exp(a + std::log(sum));
}

OMEGA_H_INLINE Real cumulative_chi_squared_density(
    Real ndofs, Real chi_squared) {
  if (chi_squared < 0.0 || ndofs < 1.0) return 0.0;
  return regularized_lower_incomplete_gamma(ndofs / 2.0, chi_squared / 2.0);
}

OMEGA_H_INLINE Real cumulative_standard_normal_density(Real x) {
  return (1.0 / 2.0) * (1.0 + std::erf(x / std::sqrt(2.0)));
}

OMEGA_H_INLINE Real cumulative_general_normal_density(
    Real mean, Real standard_deviation, Real x) {
  return cumulative_standard_normal_density((x - mean) / standard_deviation);
}

OMEGA_H_INLINE Real cumulative_weibull_density(Real shape, Real scale, Real x) {
  return 1.0 - std::exp(-std::pow(x / scale, shape));
}

}  // namespace Omega_h

#endif
