#ifndef OMEGA_H_RANDOM_HPP
#define OMEGA_H_RANDOM_HPP

#include <random123/philox.h>
#include <Omega_h_kokkos.hpp>
#include <Omega_h_few.hpp>
#include <Omega_h_scalar.hpp>

namespace Omega_h {

OMEGA_H_DEVICE Few<std::uint64_t, 2> run_philox_cbrng(Few<std::uint64_t, 2> ctr, std::uint64_t key) {
  using PhiloxType = r123::Philox4x32;
  using ctr_type = typename PhiloxType::ctr_type;
  using key_type = typename PhiloxType::key_type;
  ctr_type ctr_philox;
  key_type key_philox;
  ctr_philox[0] = static_cast<std::uint32_t>(ctr[0]);
  ctr_philox[1] = static_cast<std::uint32_t>(ctr[0] >> 32);
  ctr_philox[2] = static_cast<std::uint32_t>(ctr[1]);
  ctr_philox[3] = static_cast<std::uint32_t>(ctr[1] >> 32);
  key_philox[0] = static_cast<std::uint32_t>(key);
  key_philox[1] = static_cast<std::uint32_t>(key >> 32);
  ctr_philox = PhiloxType(ctr_philox, key_philox);
  ctr[0] = ctr_philox[1];
  ctr[0] <<= 32;
  ctr[0] |= ctr_philox[0];
  ctr[1] = ctr_philox[3];
  ctr[1] <<= 32;
  ctr[1] |= ctr_philox[2];
  return ctr;
}

OMEGA_H_INLINE Real random_uint_to_real(std::uint64_t x) {
  return static_cast<Real>(x) / static_cast<Real>(ArithTraits<std::uint64_t>::max());
}

using RandomKey = Few<std::uint64_t, 2>;
using RandomCounter = std::uint64_t;
using RandomState = Few<std::uint64_t, 2>;

OMEGA_H_DEVICE Real get_next_random_real(Few<std::uint64_t, 2> key, std::uint64_t counter,
    Few<std::uint64_t, 2>& state) {
  if (counter % 2) {
    return random_uint_to_real(state[1]);
  } else {
    Few<std::uint64_t, 2> philox_counter;
    std::uint64_t philox_key;
    philox_key = key[0];
    philox_counter[1] = key[1];
    philox_counter[0] = counter;
    state = run_philox_cbrng(key, counter);
    return random_uint_to_real(state[0]);
  }
}

}

#endif
