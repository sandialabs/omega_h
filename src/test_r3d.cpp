#include "omega_h_r3d.hpp"

#include <iostream>

int main() {
  osh::Vector<3> verts[4] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  osh::r3d::Plane<3> faces[4];
  osh::r3d::tet_faces_from_verts(faces, verts);
  for (osh::Int i = 0; i < 4; ++i) {
    std::cout << "plane #" << i << " is "
              << "(" << faces[i].n[0] << "," << faces[i].n[1] << ","
              << faces[i].n[2] << ") * " << faces[i].d << '\n';
  }
  osh::r3d::Polytope<3> a;
  osh::r3d::init_tet(&a, verts);
  std::cout << "tet nverts " << a.nverts << '\n';
  constexpr osh::Int nplanes = 1;
  osh::r3d::Plane<3> planes[nplanes] = {
    {{1,0,0},-0.5}
  };
  osh::r3d::clip(&a, planes, nplanes);
  std::cout << "after clipping, nverts = " << a.nverts << '\n';
  auto volume = osh::r3d::integrate(&a,
      osh::r3d::Polynomial<3,0>{{1}});
  std::cout << " and volume = " << volume << '\n';
}
