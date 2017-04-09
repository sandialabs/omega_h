#include "Omega_h_matrix.hpp"
#include "Omega_h_loop.hpp"

namespace Omega_h {

template <Int dim>
Reals repeat_symm(LO n, Matrix<dim, dim> symm) {
  Write<Real> symms(n * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO i) { set_symm(symms, i, symm); };
  parallel_for(n, f);
  return symms;
}

template Reals repeat_symm(LO n, Matrix<3, 3> symm);
template Reals repeat_symm(LO n, Matrix<2, 2> symm);

template <Int old_dim, Int new_dim>
Reals resize_symms_tmpl(Reals old_symms) {
  auto n = divide_no_remainder(old_symms.size(), symm_ncomps(old_dim));
  Write<Real> new_symms(n * symm_ncomps(new_dim));
  constexpr auto min_dim = min2(old_dim, new_dim);
  auto f = OMEGA_H_LAMBDA(Int i) {
    auto a = get_symm<old_dim>(old_symms, i);
    auto b = zero_matrix<new_dim, new_dim>();
    for (Int j = 0; j < min_dim; ++j) {
      for (Int k = 0; k < min_dim; ++k) {
        b[j][k] = a[j][k];
      }
    }
    set_symm(new_symms, i, b);
  };
  parallel_for(n, f);
  return new_symms;
}

Reals resize_symms(Reals old_symms, Int old_dim, Int new_dim) {
  if (old_dim == new_dim) return old_symms;
  if (old_dim == 2 && new_dim == 3) return resize_symms_tmpl<2, 3>(old_symms);
  if (old_dim == 3 && new_dim == 2) return resize_symms_tmpl<3, 2>(old_symms);
  OMEGA_H_NORETURN(Reals());
}

}  // end namespace Omega_h
