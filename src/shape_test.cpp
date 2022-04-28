#include "Omega_h_shape.hpp"
#include <iostream>


void check_verts(const Omega_h::Matrix<2,3> & coords) {
  for(int i=0; i<coords.size(); ++i ) {
    auto xi = Omega_h::barycentric_from_global<2,2>(coords[i], coords);
    Omega_h::Vector<3> hand_xi{0,0,0};
    hand_xi[i] = 1;
    printf("[%f,%f,%f] == [%f,%f,%f]\n",xi[0],xi[1],xi[2],hand_xi[0],hand_xi[1],hand_xi[2]);
    OMEGA_H_CHECK(Omega_h::are_close(xi,hand_xi));
  }
}

void check_point(const Omega_h::Matrix<2,3> & coords) {
    Omega_h::Vector<2> point{0.5,0.5};
    auto xi = Omega_h::barycentric_from_global<2,2>(point, coords);
    Omega_h::Vector<3> hand_xi{0.25,0.25,0.5};
    printf("[%f,%f,%f] == [%f,%f,%f]\n",xi[0],xi[1],xi[2],hand_xi[0],hand_xi[1],hand_xi[2]);
    OMEGA_H_CHECK(are_close(xi,hand_xi));

}


int main(int argc, char** argv) {
  Omega_h::Matrix<2,3> coords{{0.0},{1,0},{0.5,1}} ;
  check_verts(coords);
  check_point(coords);

  return 0;
}
