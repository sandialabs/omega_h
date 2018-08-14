#include <cmath>
#include <iostream>
#include <numeric>

#include "Omega_h_file.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_mesh.hpp"

/*!
 *  \brief This program computes the gradient of the same function
 *  as the one defined in \e field_on_square example.
 *  \author Omar Awile
 *  \author Samuel Melchior
*/
int main(int argc, char** argv) {
  // initialization
  auto lib = Omega_h::Library(&argc, &argv);
  const auto world = lib.world();
  auto mesh = Omega_h::gmsh::read("square.msh", world);
  const auto tris2verts = mesh.ask_elem_verts();
  const auto coords = mesh.coords();

  // synchronized parallel print of the number of vertices on each proc
  const auto rank = mesh.comm()->rank();
  const auto nranks = mesh.comm()->size();
  for (int turn = 0; turn < nranks; turn++) {
    if (rank == turn) {
      std::cout << "[" << rank << "] number of vertices " << mesh.nverts()
                << std::endl;
    }
    world->barrier();
  }

  // same quadratic function as in fieldOnSquare client
  const Omega_h::Write<Omega_h::Real> u_w(mesh.nverts());
  const auto init_u = OMEGA_H_LAMBDA(Omega_h::LO r) {
    auto x_r = Omega_h::get_vector<2>(coords, r);
    u_w[r] = x_r[1] - std::pow(2 * x_r[0] - 1, 2);
  };
  Omega_h::parallel_for(mesh.nverts(), init_u);
  const Omega_h::Reals u_r(u_w);
  mesh.add_tag(Omega_h::VERT, "u", 1, u_r);

  // the size of this array must be twice the number of elements
  const Omega_h::Write<Omega_h::Real> gradu_w(mesh.nelems() * 2);
  const auto calc_gradient = OMEGA_H_LAMBDA(Omega_h::LO j) {

    const auto tri_j2verts = Omega_h::gather_verts<3>(tris2verts, j);

    // coords (x, y) and value of u at the 3 vertices in the jth triangle
    const auto tri_j2x = Omega_h::gather_vectors<3, 2>(coords, tri_j2verts);
    const auto tri_j2u = Omega_h::gather_scalars<3>(u_r, tri_j2verts);

    Omega_h::Matrix<2, 2> M;
    Omega_h::Vector<2> b;
    for (auto ir = 0; ir < 2; ++ir) {
      for (auto ic = 0; ic < 2; ++ic) {
        M[ic][ir] = tri_j2x[ir + 1][ic] - tri_j2x[0][ic];
      }
      b[ir] = tri_j2u[ir + 1] - tri_j2u[0];
    }
    const auto tri_j2grad = invert(M) * b;
    gradu_w[2 * j] = tri_j2grad[0];
    gradu_w[2 * j + 1] = tri_j2grad[1];
  };
  Omega_h::parallel_for(mesh.nelems(), calc_gradient);
  // add a 2d vector field with values stored at each triangle
  mesh.add_tag(mesh.dim(), "grad_u", 2, Omega_h::Reals(gradu_w));

  Omega_h::vtk::write_parallel("gradientField", &mesh);
}
