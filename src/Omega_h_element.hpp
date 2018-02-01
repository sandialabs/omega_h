#ifndef OMEGA_H_ELEMENT_HPP
#define OMEGA_H_ELEMENT_HPP

#include <Omega_h_c.h>
#include <Omega_h_simplex.hpp>
#include <Omega_h_hypercube.hpp>

namespace Omega_h {

OMEGA_H_INLINE Int element_down_template(Omega_h_Family family,
    Int elem_dim, Int bdry_dim, Int which_bdry, Int which_vert) {
  return (family == OMEGA_H_SIMPLEX ?
          simplex_down_template(elem_dim, bdry_dim, which_bdry, which_vert) :
          hypercube_down_template(elem_dim, bdry_dim, which_bdry, which_vert));
}

OMEGA_H_INLINE Int element_degree(Omega_h_Family family, Int from_dim, Int to_dim) {
  return (family == OMEGA_H_SIMPLEX ?
          simplex_degree(from_dim, to_dim) :
          hypercube_degree(from_dim, to_dim));
}

OMEGA_H_INLINE TemplateUp element_up_template(
    Omega_h_Family family, Int elem_dim, Int bdry_dim, Int which_bdry, Int which_up) {
  return (family == OMEGA_H_SIMPLEX ?
          simplex_up_template(elem_dim, bdry_dim, which_bdry, which_up) :
          hypercube_up_template(elem_dim, bdry_dim, which_bdry, which_up));
}

}

#endif
