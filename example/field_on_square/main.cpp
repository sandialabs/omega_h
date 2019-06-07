#include <cmath>

#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"

/*!
 *  \brief Minimalist Omega_h client
 *  \details Load a square mesh, looping on whole its vertices to
 *  assign the value y-(2x-1)² from the coordinates (x,y) and output it in
 *  a vtu file.
 *  \author Samuel Melchior
*/
int main(int argc, char** argv) {
  // calls MPI_init(&argc, &argv) and Kokkos::initialize(argc, argv)
  auto lib = Omega_h::Library(&argc, &argv);
  // encapsulates many MPI functions, e.g., world.barrier()
  const auto world = lib.world();
  // create a distributed mesh object, calls mesh.balance() and then returns it
  auto mesh = Omega_h::gmsh::read("square.msh", world);
  // return type Omega_h::Reals derived from Omega_h::Read<Omega_h::Real>
  const auto coords = mesh.coords();
  // array of type Omega_h::Write must be used when modified as below
  const Omega_h::Write<Omega_h::Real> u_w(mesh.nverts());
  // Omega_h::LO is a 32 bits int for local indexes on each proc
  const auto initialize_u = OMEGA_H_LAMBDA(Omega_h::LO r) {
    // get_vector<2> abstracts the storage convention inside coords array
    const auto x_r = Omega_h::get_vector<2>(coords, r);
    // quadratic 2d function between -1 and 1 on the square [0,1] x [0,1]
    u_w[r] = x_r[1] - std::pow(2 * x_r[0] - 1, 2);
  };
  // encapsulates Kokkos::parallel_for
  Omega_h::parallel_for(mesh.nverts(), initialize_u);
  // create a read_only array u_r from the initialized field "written" in u_w
  Omega_h::Reals u_r(u_w);
  // finally enforcing that u_r correspond to a scalar field "u" at VERTices
  mesh.add_tag(Omega_h::VERT, "u", 1, u_r);
  // export the results to visualize in paraview
  Omega_h::vtk::write_vtu("field_on_square.vtu", &mesh, mesh.dim());
}
