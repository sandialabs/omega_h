#include "algebra.hpp"
#include "access.hpp"
#include "loop.hpp"

namespace Omega_h {

template <Int dim>
static Reals get_vector_norms_tmpl(Reals vs) {
  CHECK(vs.size() % dim == 0);
  auto n = vs.size() / dim;
  auto out = Write<Real>(n);
  auto f = LAMBDA(LO i) { out[i] = norm(get_vector<dim>(vs, i)); };
  parallel_for(n, f);
  return out;
}

Reals get_vector_norms(Reals vs, Int dim) {
  switch (dim) {
    case 3:
      return get_vector_norms_tmpl<3>(vs);
    case 2:
      return get_vector_norms_tmpl<2>(vs);
  }
  NORETURN(Reals());
}

template <Int dim>
static Reals normalize_vectors_tmpl(Reals vs) {
  CHECK(vs.size() % dim == 0);
  auto n = vs.size() / dim;
  auto out = Write<Real>(vs.size());
  auto f = LAMBDA(LO i) {
    set_vector<dim>(out, i, normalize(get_vector<dim>(vs, i)));
  };
  parallel_for(n, f);
  return out;
}

Reals normalize_vectors(Reals vs, Int dim) {
  switch (dim) {
    case 3:
      return normalize_vectors_tmpl<3>(vs);
    case 2:
      return normalize_vectors_tmpl<2>(vs);
  }
  NORETURN(Reals());
}

Reals interpolate_between(Reals a, Reals b, Real t) {
  CHECK(a.size() == b.size());
  auto n = a.size();
  auto out = Write<Real>(n);
  auto f = LAMBDA(LO i) { out[i] = a[i] * (1.0 - t) + b[i] * t; };
  parallel_for(n, f);
  return out;
}
}  // namespace Omega_h
