#include "Omega_h_eigen.hpp"

#include "Omega_h_loop.hpp"

namespace Omega_h {

template <Int dim>
Reals get_max_eigenvalues_dim(Reals symms) {
  auto n = divide_no_remainder(symms.size() / symm_ncomps(dim));
  auto out = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto a = get_symm<dim>(symms, i);
    auto ews = get_eigenvalues(a);
    auto best_ew = fabs(ews.values[0]);
    for (Int j = 1; j < ews.n; ++j) {
      auto cand = fabs(ews.values[j]);
      if (cand > best_ew) best_ew = cand;
    }
    out[i] = best_ew;
  };
  parallel_for(n, f);
  return out;
}

Reals get_max_eigenvalues(Int dim, Reals symms) {
  if (dim == 3) return get_max_eigenvalues_dim<3>(symms);
  if (dim == 2) return get_max_eigenvalues_dim<2>(symms);
  OMEGA_H_NORETURN(Reals());
}

}  // end namespace Omega_h
