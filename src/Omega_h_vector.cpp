#include "Omega_h_vector.hpp"
#include "Omega_h_fail.hpp"
#include "Omega_h_for.hpp"

namespace Omega_h {

template <Int dim>
Reals get_vector_norms_tmpl(Reals vs) {
  auto n = divide_no_remainder(vs.size(), dim);
  auto out = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) { out[i] = norm(get_vector<dim>(vs, i)); };
  parallel_for(n, f, "get_vector_norms");
  return out;
}

Reals get_vector_norms(Reals vs, Int dim) {
  switch (dim) {
    case 3:
      return get_vector_norms_tmpl<3>(vs);
    case 2:
      return get_vector_norms_tmpl<2>(vs);
  }
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals normalize_vectors_tmpl(Reals vs) {
  auto n = divide_no_remainder(vs.size(), dim);
  auto out = Write<Real>(vs.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    set_vector<dim>(out, i, normalize(get_vector<dim>(vs, i)));
  };
  parallel_for(n, f, "normalize_vectors");
  return out;
}

Reals normalize_vectors(Reals vs, Int dim) {
  switch (dim) {
    case 3:
      return normalize_vectors_tmpl<3>(vs);
    case 2:
      return normalize_vectors_tmpl<2>(vs);
  }
  OMEGA_H_NORETURN(Reals());
}

template <Int dim>
Reals dot_vectors_dim(Reals a, Reals b) {
  OMEGA_H_CHECK(a.size() == b.size());
  auto n = divide_no_remainder(a.size(), dim);
  auto out = Write<Real>(n);
  auto f = OMEGA_H_LAMBDA(LO i) {
    out[i] = get_vector<dim>(a, i) * get_vector<dim>(b, i);
  };
  parallel_for(n, f, "dot_vectors");
  return out;
}

Reals dot_vectors(Reals a, Reals b, Int dim) {
  if (dim == 3) return dot_vectors_dim<3>(a, b);
  if (dim == 2) return dot_vectors_dim<2>(a, b);
  if (dim == 1) return dot_vectors_dim<1>(a, b);
  OMEGA_H_NORETURN(Reals());
}

Reals resize_vectors(Reals old_vectors, Int old_dim, Int new_dim) {
  if (old_dim == new_dim) return old_vectors;
  auto nv = divide_no_remainder(old_vectors.size(), old_dim);
  Write<Real> new_vectors(nv * new_dim);
  auto min_dim = min2(old_dim, new_dim);
  auto f = OMEGA_H_LAMBDA(Int i) {
    for (Int j = 0; j < min_dim; ++j) {
      new_vectors[i * new_dim + j] = old_vectors[i * old_dim + j];
    }
    for (Int j = min_dim; j < new_dim; ++j) {
      new_vectors[i * new_dim + j] = 0.0;
    }
  };
  parallel_for(nv, f, "resize_vectors");
  return new_vectors;
}

template <Int dim>
Reals repeat_vector(LO n, Vector<dim> v) {
  Write<Real> vs(n * dim);
  auto f = OMEGA_H_LAMBDA(LO i) { set_vector(vs, i, v); };
  parallel_for(n, f, "repeat_vector");
  return vs;
}

template Reals repeat_vector(LO n, Vector<1> v);
template Reals repeat_vector(LO n, Vector<2> v);
template Reals repeat_vector(LO n, Vector<3> v);
template Reals repeat_vector(LO n, Vector<4> v);
template Reals repeat_vector(LO n, Vector<6> v);
template Reals repeat_vector(LO n, Vector<9> v);

}  // end namespace Omega_h
