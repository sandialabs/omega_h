#include <cmath>
#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_loop.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mark.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_shape.hpp"

/*!
 *  \brief Computes finite element discretization of the 2nd Poisson equation
 *  \details This program computes the solution of the finite element
 *  discretization of the 2nd Poisson equation with constant 
 *  right hand side = 3 and then evaluates the local residual at each vertex
 *  of the corresponding finite volume discretization.
 *  \remark the iterative solution is a Richardson iteration
 *  with a <b>Jacobi preconditioner</b>.
 *  \author Samuel Melchior
 */
int main(int argc, char** argv) {
  // initialization
  auto lib = Omega_h::Library(&argc, &argv);
  const auto world = lib.world();
  auto mesh = Omega_h::gmsh::read("square.msh", world);

  // very important to call partitioning with ghost layer first 
  mesh.set_parting(OMEGA_H_GHOSTED);
  const auto tris2verts = mesh.ask_elem_verts();
  const auto coords = mesh.coords();
  const auto area = measure_elements_real(&mesh);

  // computing global and local information on number and area of triangles
  const auto dim = mesh.dim();
  const auto correctGlobalArea =
      get_sum(mesh.comm(), mesh.owned_array(dim, area, 1));
  const auto correctOwnedArea = get_sum(mesh.owned_array(dim, area, 1));
  const auto wrongTotalArea = get_sum(mesh.comm(), area);
  const auto wrongOwnedArea = get_sum(area);
  const auto nGlobalTriangles = mesh.nglobal_ents(dim);
  const auto nOwnedTriangles = get_sum(mesh.owned(dim));

  // synchronized parallel print including information specific to ghost elements
  const auto rank = mesh.comm()->rank();
  if (rank == 0) {
    std::cout << "Globally, " << nGlobalTriangles << " cover a total area of "
              << correctGlobalArea << " and not " << wrongTotalArea
              << std::endl;
  }
  const auto nranks = mesh.comm()->size();
  for (auto turn = 0; turn < nranks; turn++) {
    if (rank == turn)
      std::cout << "On rank " << rank << ", " << nOwnedTriangles
                << " cover a local area of " << correctOwnedArea << " ["
                << wrongOwnedArea - correctOwnedArea << " covered by "
                << mesh.nelems() - nOwnedTriangles << " ghost triangles]"
                << std::endl;
    world->barrier();
  }

  // initialize the array for the candidate solution with value 0
  const Omega_h::Write<Omega_h::Real> sol_w(mesh.nverts(), 0);
  // read-only arrays for the current residual and candidate solution
  Omega_h::Reals res, sol;
  // on the current mesh, convergence to machine precision in < 1000 iterations
  for (auto iter = 0; iter < 1000; iter++) {
    // initialization to 0 before assembling the contribution of each triangle
    const Omega_h::Write<Omega_h::Real> res_w(mesh.nverts(), 0);
    const Omega_h::Write<Omega_h::Real> diag_w(mesh.nverts(), 0);

    sol = mesh.sync_array(Omega_h::VERT, Omega_h::Reals(sol_w), 1);

    const auto computeResidual = OMEGA_H_LAMBDA(Omega_h::LO j) {
      // topology and geometry of triangle with local index j
      const auto tri_j2verts = Omega_h::gather_verts<3>(tris2verts, j);
      const auto tri_j2x = Omega_h::gather_vectors<3, 2>(coords, tri_j2verts);

      // M[ic][ir] = tri_j2x[ir + 1][ic] - tri_j2x[0][ic] in gradient/main.cpp
      const auto M = Omega_h::simplex_basis<2, 2>(tri_j2x);
      // value of the current solution at the vertices in tri_j2verts
      const auto tri_j2sol = Omega_h::gather_scalars<3>(sol, tri_j2verts);

      // projection matrix and solution vector such that u*P = b in gradient
      const auto P = Omega_h::Matrix<3, 2>({{-1, 1, 0}, {-1, 0, 1}});
      Omega_h::Vector<3> u;
      for (auto ir = 0; ir < 3; ++ir) {
        u[ir] = tri_j2sol[ir];
      }

      // u*grad_phi = grad_u is the gradient of the finite element interpolation
      const auto grad_phi = P * invert(M);
      // grad_phi[ic][ir] = dphi_r/dx_c
      const auto localStiffness = -grad_phi * transpose(grad_phi) * area[j];
      // integration by parts of the Laplacian yields a minus sign
      const auto local_lhs = localStiffness * u;
      for (auto ir = 0; ir < 3; ++ir) {
        // continuous rhs = -3 because \int phi = 1/6 on parent element
        res_w[tri_j2verts[ir]] -= local_lhs[ir] + area[j];
        diag_w[tri_j2verts[ir]] += localStiffness[ir][ir];
      }
    };
    Omega_h::parallel_for(mesh.nelems(), computeResidual);

    // loop on vertices that are on the four edges of the square geometry
    for (std::size_t i = 0; i < 4; ++i) {
      // return type Omega_h::Read<Omega_h::I8> to store booleans
      const auto rs_are_brs =
          mark_class_closure(&mesh, Omega_h::VERT, Omega_h::EDGE, i + 1);
      // br+1 is a boundary index between 1 and the number of boundary vertices
      const auto br2r = collect_marked(rs_are_brs);

      // residual = 0 is sufficient to impose a Dirichlet boundary condition
      const auto applyHomogeneousDirichletBC = OMEGA_H_LAMBDA(Omega_h::LO br) {
        res_w[br2r[br]] = 0;
      };
      Omega_h::parallel_for(br2r.size(), applyHomogeneousDirichletBC);
    }
    // ghosting ensures the residual is correctly computed on owned vertices
    res = mesh.sync_array(Omega_h::VERT, Omega_h::Reals(res_w), 1);

    // monitoring the convergence of the residual
    const auto localSquaredRes = multiply_each(res, res);
    const auto maxRes =
        sqrt(get_max(mesh.comm(), Omega_h::Reals(localSquaredRes)));
    if (rank == 0 && iter % 50 == 0) {
      std::cout << "[" << iter << "]: " << maxRes << std::endl;
    }

    // Richardson iteration with a Jacobi preconditioner
    const auto updateSolution = OMEGA_H_LAMBDA(Omega_h::LO r) {
      sol_w[r] += 0.5 * res_w[r] / diag_w[r];
    };
    Omega_h::parallel_for(mesh.nverts(), updateSolution);
  }

  // export residual and solution as scalar fields at the vertices
  mesh.add_tag(Omega_h::VERT, "res", 1, res);
  mesh.add_tag(Omega_h::VERT, "sol", 1, sol);
  Omega_h::vtk::write_parallel("Poisson", &mesh);
}
