#include <Omega_h_for.hpp>
#include <Omega_h_random.hpp>
#include <Omega_h_random_inline.hpp>

namespace Omega_h {

Reals unit_uniform_random_reals_from_globals(
    GOs const globals, I64 const seed, I64 const counter) {
  auto const n = globals.size();
  auto const out = Write<Real>(n);
  auto functor = OMEGA_H_LAMBDA(LO const i) {
    auto const global = globals[i];
    UnitUniformDistribution distrib(seed, global, counter);
    out[i] = distrib();
  };
  parallel_for(n, std::move(functor));
  return out;
}

}  // namespace Omega_h
