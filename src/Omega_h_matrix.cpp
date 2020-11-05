#include "Omega_h_matrix.hpp"
#include "Omega_h_for.hpp"

namespace Omega_h {

template <Int dim>
Reals repeat_symm(LO const n, Tensor<dim> const symm) {
  return repeat_vector(n, symm2vector(symm));
}

template Reals repeat_symm(LO const n, Tensor<3> const symm);
template Reals repeat_symm(LO const n, Tensor<2> const symm);
template Reals repeat_symm(LO const n, Tensor<1> const symm);

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
  parallel_for(n, f, "resize_symms");
  return new_symms;
}

Reals resize_symms(Reals old_symms, Int old_dim, Int new_dim) {
  if (old_dim == new_dim) return old_symms;
  if (old_dim == 2 && new_dim == 3) return resize_symms_tmpl<2, 3>(old_symms);
  if (old_dim == 3 && new_dim == 2) return resize_symms_tmpl<3, 2>(old_symms);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals repeat_matrix(LO const n, Tensor<dim> const m) {
  return repeat_vector(n, matrix2vector(m));
}

template Reals repeat_matrix(LO const n, Tensor<3> const m);
template Reals repeat_matrix(LO const n, Tensor<2> const m);
template Reals repeat_matrix(LO const n, Tensor<1> const m);

template <Int dim>
Reals matrices_times_vectors_dim(Reals ms, Reals vs) {
  auto n = divide_no_remainder(vs.size(), dim);
  OMEGA_H_CHECK(ms.size() == n * matrix_ncomps(dim, dim));
  auto out = Write<Real>(n * dim);
  auto f = OMEGA_H_LAMBDA(LO i) {
    set_vector(out, i, get_matrix<dim, dim>(ms, i) * get_vector<dim>(vs, i));
  };
  parallel_for(n, f, "matrices_times_vectors");
  return out;
}

Reals matrices_times_vectors(Reals ms, Reals vs, Int dim) {
  if (dim == 3) return matrices_times_vectors_dim<3>(ms, vs);
  if (dim == 2) return matrices_times_vectors_dim<2>(ms, vs);
  if (dim == 1) return matrices_times_vectors_dim<1>(ms, vs);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals matrices_times_matrices_dim(Reals a, Reals b) {
  auto n = divide_no_remainder(a.size(), matrix_ncomps(dim, dim));
  OMEGA_H_CHECK(b.size() == n * matrix_ncomps(dim, dim));
  auto out = Write<Real>(n * matrix_ncomps(dim, dim));
  auto f = OMEGA_H_LAMBDA(LO i) {
    set_matrix(out, i, get_matrix<dim, dim>(a, i) * get_matrix<dim, dim>(b, i));
  };
  parallel_for(n, f, "matrices_times_matrices");
  return out;
}

Reals matrices_times_matrices(Reals a, Reals b, Int dim) {
  if (dim == 3) return matrices_times_matrices_dim<3>(a, b);
  if (dim == 2) return matrices_times_matrices_dim<2>(a, b);
  if (dim == 1) return matrices_times_matrices_dim<1>(a, b);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals symms_inria2osh_dim(Reals symms) {
  auto n = divide_no_remainder(symms.size(), symm_ncomps(dim));
  Write<Real> out(symms.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto iv = get_vector<symm_ncomps(dim)>(symms, i);
    auto is = vector2symm_inria(iv);
    auto ov = symm2vector(is);
    set_vector(out, i, ov);
  };
  parallel_for(n, f);
  return out;
}

Reals symms_inria2osh(Int dim, Reals symms) {
  if (dim == 3) return symms_inria2osh_dim<3>(symms);
  if (dim == 2) return symms_inria2osh_dim<2>(symms);
  if (dim == 1) return symms_inria2osh_dim<1>(symms);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals symms_osh2inria_dim(Reals symms) {
  auto n = divide_no_remainder(symms.size(), symm_ncomps(dim));
  Write<Real> out(symms.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    auto iv = get_vector<symm_ncomps(dim)>(symms, i);
    auto is = vector2symm(iv);
    auto ov = symm2vector_inria(is);
    set_vector(out, i, ov);
  };
  parallel_for(n, f);
  return out;
}

Reals symms_osh2inria(Int dim, Reals symms) {
  if (dim == 3) return symms_osh2inria_dim<3>(symms);
  if (dim == 2) return symms_osh2inria_dim<2>(symms);
  if (dim == 1) return symms_osh2inria_dim<1>(symms);
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals matrices_to_symms_dim(Reals const matrices) {
  constexpr auto ncomps_in = square(dim);
  constexpr auto ncomps_out = symm_ncomps(dim);
  auto const n = divide_no_remainder(matrices.size(), ncomps_in);
  auto const out = Write<Real>(n * ncomps_out);
  auto functor = OMEGA_H_LAMBDA(LO const i) {
    set_symm(out, i, get_matrix<dim, dim>(matrices, i));
  };
  parallel_for(n, std::move(functor));
  return out;
}

Reals matrices_to_symms(Reals const matrices, Int const dim) {
  if (dim == 3) return matrices_to_symms_dim<3>(matrices);
  if (dim == 2) return matrices_to_symms_dim<2>(matrices);
  if (dim == 1) return matrices_to_symms_dim<1>(matrices);
  OMEGA_H_NORETURN(Reals());
}

}  // end namespace Omega_h
