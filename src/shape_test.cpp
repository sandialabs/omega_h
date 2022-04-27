#include "Omega_h_shape.hpp"
#include <iostream>


auto get_inverse_basis(const Omega_h::Matrix<2,3> & coords) { 
  auto basis = Omega_h::simplex_basis<2,2>(coords);
  Omega_h::Matrix<2,2> hand_basis {{coords(0,1)-coords(0,0),coords(1,1)-coords(1,0)},
                                   {coords(0,2)-coords(0,0),coords(1,2)-coords(1,0)}};
  OMEGA_H_CHECK(are_close(basis, hand_basis));
  auto inverse_basis = Omega_h::invert(basis);
  return inverse_basis;
}


void check_verts(const Omega_h::Matrix<2,3> & coords, const Omega_h::Matrix<2,2> & inverse_basis) {
  for(int i=0; i<coords.size(); ++i ) {
    auto lambda = inverse_basis*(coords[i] - coords[0]);
    auto xi = Omega_h::form_barycentric(lambda);
    Omega_h::Vector<3> hand_xi{0,0,0};
    hand_xi[i] = 1;
    OMEGA_H_CHECK(Omega_h::are_close(xi,hand_xi));
  }
}

void check_point(const Omega_h::Matrix<2,3> & coords, const Omega_h::Matrix<2,2> & inverse_basis) {
    Omega_h::Vector<2> point{0.5,0.5};
    auto lambda = inverse_basis*(point - coords[0]);
    auto xi = Omega_h::form_barycentric(lambda);

    Omega_h::Vector<3> hand_xi{0.25,0.25,0.5};
    OMEGA_H_CHECK(are_close(xi,hand_xi));

}


int main(int argc, char** argv) {
  Omega_h::Matrix<2,3> coords{{0.0},{1,0},{0.5,1}} ;
  auto inverse_basis = get_inverse_basis(coords);
  check_verts(coords, inverse_basis);
  check_point(coords,inverse_basis);

  return 0;
}
