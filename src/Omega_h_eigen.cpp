#include "Omega_h_eigen.hpp"

#include "Omega_h_for.hpp"

namespace Omega_h {

template <Int n_max, typename T>
OMEGA_H_INLINE T maximum_magnitude(Few<T, n_max> x, Int n) {
  auto out = x[0];
  auto max_mag = std::abs(x[0]);
  for (Int i = 1; i < n; ++i) {
    auto mag = std::abs(x[i]);
    if (mag > max_mag) {
      max_mag = mag;
      out = x[i];
    }
  }
  return out;
}

template <Int dim>
Reals get_max_eigenvalues_dim(Reals symms) {
  auto n = divide_no_remainder(symms.size(), symm_ncomps(dim));
  auto out = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto a = get_symm<dim>(symms, i);
    auto ews = get_eigenvalues(a);
    auto max_ew = maximum_magnitude(ews.values, ews.n);
    out[i] = max_ew;
  };
  parallel_for(n, f, "get_max_eigenvalues");
  return out;
}

Reals get_max_eigenvalues(Int dim, Reals symms) {
  if (dim == 3) return get_max_eigenvalues_dim<3>(symms);
  if (dim == 2) return get_max_eigenvalues_dim<2>(symms);
  if (dim == 1) return symms;
  OMEGA_H_NORETURN(Reals());
}

}  // end namespace Omega_h
