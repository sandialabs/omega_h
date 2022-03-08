#include "Omega_h_build.hpp"
#include "Omega_h_library.hpp"
#include <Omega_h_file.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_remotes.hpp>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  assert(world->size() == 2);

  const int dim = 2;
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, 1., 1., 0., 2, 2, 0);
  assert(mesh.nglobal_ents(dim) == 8);
  assert(mesh.dim() == dim);
  Omega_h::vtk::write_parallel("box_before.vtk", &mesh, dim);
  const int rank = world->rank();

  //rank 0 has elms 1,2,4,5
  //rank 1 has elms 0,3,6,7
  //move elm with gid 7 from rank 1 to rank 0
  const int elms = (!rank) ? 5 : 3;
  if(!rank) {
    Write<I32> ranks = {0,0,0,0,1};
    Write<LO> idxs = {0,1,2,3,3};
    auto owners = Remotes(ranks, idxs);
    mesh.migrate(owners);
  } else {
    Write<I32> ranks = {1,1,1};
    Write<LO> idxs = {0,1,2};
    auto owners = Remotes(ranks, idxs);
    mesh.migrate(owners);
  }
  if(!rank) {
    OMEGA_H_CHECK(mesh.nelems()==5);
  } else {
    OMEGA_H_CHECK(mesh.nelems()==3);
  }
  Omega_h::vtk::write_parallel("box_after.vtk", &mesh, mesh.dim());
  return 0;
}

