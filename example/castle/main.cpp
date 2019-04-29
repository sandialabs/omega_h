#include <cmath>
#include <iostream>
#include <numeric>

#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_mesh.hpp"

/*!
 *  \brief Sandbox castle
 *  \details This program is intended to be run with `mpirun -n 2 ./castle`.
 *  The local and global indices of the vertices of each triangle
 *  in a very coarse mesh are written to standard output.
 *  Then some iterations over the triangles and the vertices are performed
 *  to investigate the usage of ghost elements depending
 *  on the definition of the `NOGHOSTS` preprocessor macro.
 *  \author Omar Awile
 *  \author Samuel Melchior
 */

/*
// preprocessor variable that prevents the creation of ghost elements
#define NOGHOSTS
*/

int main(int argc, char** argv) {
  // initialization of MPI, Kokkos and the distributed mesh
  auto lib = Omega_h::Library(&argc, &argv);
  const auto world = lib.world();
  auto mesh = Omega_h::gmsh::read("square.msh", world);

  // return the topology from each triangle to its 3 vertices as Omega_h::LOs
  auto tris2verts = mesh.ask_elem_verts();

  // get MPI rank of this proc and total number of MPI procs
  const auto rank = mesh.comm()->rank();
  const auto nranks = mesh.comm()->size();

  // return type GOs is an array of 64bits long int for global indices
  const auto globals_v = mesh.globals(Omega_h::VERT);
  const auto globals_e = mesh.globals(Omega_h::FACE);

  // synchronized parallel print to avoid interlacing between procs
  for (int turn = 0; turn < nranks; turn++) {
    // it is only the turn of one processor to print its information
    if (rank == turn) {
      std::cout << "[" << rank << "] number of vertices " << mesh.nverts()
                << std::endl;
      for (int local_index = 0; local_index < mesh.nverts(); ++local_index) {
        // global_v maps local_index on this proc to global index for vertices
        std::cout << "[" << rank << "]: "
                  << "Local vertex number " << local_index
                  << " has global index " << globals_v[local_index]
                  << std::endl;
      }
      for (int j = 0; j < mesh.nelems(); ++j) {
        // return type Omega_h::Few<Omega_h::LO, 3> encapsulates a static array
        const auto tri_j2verts = Omega_h::gather_verts<3>(tris2verts, j);
        // tri_j2verts gives local index on this proc of each triangle vertex
        std::cout << "[" << rank << "]: "
                  << "local  " << j << " " << tri_j2verts[0] << " "
                  << tri_j2verts[1] << " " << tri_j2verts[2] << std::endl;
        std::cout << "[" << rank << "]: "
                  << "global " << globals_e[j] << " "
                  << globals_v[tri_j2verts[0]] << " "
                  << globals_v[tri_j2verts[1]] << " "
                  << globals_v[tri_j2verts[2]] << std::endl;
      }
    }
    world->barrier();
  }

#ifndef NOGHOSTS
  // adding one layer of ghost elements for correct communication between procs
  mesh.set_parting(OMEGA_H_GHOSTED);
  // call this after creating ghostlayers to update this data structure
  tris2verts = mesh.ask_elem_verts();
#endif

  // create a writeable data array and initialize it to the current rank
  const auto vert_rank =
      Omega_h::Reals(Omega_h::Write<Omega_h::Real>(mesh.nverts(), rank));
  mesh.add_tag(Omega_h::VERT, "rank", 1,
      vert_rank);  // register the data array in the mesh
  auto request = mesh.isync_array(Omega_h::VERT, vert_rank, 1);  // update ghosts
  auto finished = request.completed(); // ask for synchronization completion
  if (!finished) {
    // do work
  }
  const auto vert_rank_s = request.get(); // wait and get result
  mesh.add_tag(
      Omega_h::VERT, "synced_rank", 1, vert_rank_s);  // synchronized rank

  // create a copy of type Omega_h::Write<Omega_h::Real> from a Reals array
  const auto vrank_srw = Omega_h::deep_copy(vert_rank_s);
  // the value at each vertex is increased by 20% of each neighboring vertex
  const auto elements2vertices = OMEGA_H_LAMBDA(Omega_h::LO j) {
    const auto tri_j2verts = Omega_h::gather_verts<3>(tris2verts, j);
    double sumVerts = 0.;
    for (auto i = 0; i < 3; ++i) {
      sumVerts += vert_rank_s[tri_j2verts[i]];
    }
    for (auto i = 0; i < 3; ++i) {
      vrank_srw[tri_j2verts[i]] += 0.1 * sumVerts;
    }
  };
  Omega_h::parallel_for(mesh.nelems(), elements2vertices);
  mesh.add_tag(Omega_h::VERT, "weighted_sum", 1, Omega_h::Reals(vrank_srw));
  const auto vrank_srw_sync =
      mesh.sync_array(Omega_h::VERT, Omega_h::Reals(vrank_srw), 1);
  mesh.add_tag(Omega_h::VERT, "synced_weight_sum", 1, vrank_srw_sync);

// export results into pieces.pvtu for parallel visualization in paraview
#ifndef NOGHOSTS
  Omega_h::vtk::write_parallel("ghosted", &mesh);
#else
  Omega_h::vtk::write_parallel("noGhost", &mesh);
#endif
}
