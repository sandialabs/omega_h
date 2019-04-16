#include <Omega_h_amr.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <iostream>

namespace Omega_h {

static void test_2D_case1(CommPtr comm) {
  auto mesh = build_box(comm, OMEGA_H_HYPERCUBE,
      2, 1., 0., 2, 1, 0); 
  vtk::FullWriter writer("out_amr_mpi_2D", &mesh);
  writer.write();
  {
    Write<Byte> marks(mesh.nelems(), 0);
    if (comm->rank() == 1) marks.set(0, 1);
    const auto xfer_opts = TransferOpts();
    amr::refine(&mesh, marks, xfer_opts);
  }
  if (comm->rank() == 0) {
    std::cout << "nelems on rank 0 after first adapt: " << mesh.nelems() << '\n';
    std::cout << "nedges on rank 0 after first adapt: " << mesh.nedges() << '\n';
    std::cout << "nverts on rank 0 after first adapt: " << mesh.nverts() << '\n';
  }
  writer.write();
  {
    Write<Byte> marks(mesh.nelems(), 0); 
    if (comm->rank() == 1) marks.set(2, 1);
    const auto xfer_opts = TransferOpts();
    amr::refine(&mesh, marks, xfer_opts);
  }
  // why did the two child edges and one child vert on rank 0
  // after the second adapt when no elements on the partition
  // boundary were modified?
  if (comm->rank() == 0) {
    std::cout << "nelems on rank 0 after second adapt: " << mesh.nelems() << '\n';
    std::cout << "nedges on rank 0 after second adapt: " << mesh.nedges() << '\n';
    std::cout << "nverts on rank 0 after second adapt: " << mesh.nverts() << '\n';
  }
  writer.write();
}

static void test_2D_case2(CommPtr comm) {
  auto mesh = build_box(comm, OMEGA_H_HYPERCUBE,
      2, 1., 0., 2, 1, 0); 
  { // refine element on rank 0
    Write<Byte> marks(mesh.nelems(), 0);
    if (comm->rank() == 0) marks.set(0, 1);
    const auto xfer_opts = TransferOpts();
    amr::refine(&mesh, marks, xfer_opts);
  }
  { // refine element on rank 1 (dies in get_amr_topology)
    Write<Byte> marks(mesh.nelems(), 0);
    if (comm->rank() == 1) marks.set(0, 1);
    const auto xfer_opts = TransferOpts();
    amr::refine(&mesh, marks, xfer_opts);
  }
}

}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  OMEGA_H_CHECK(comm->size() == 2); 
  Omega_h::test_2D_case1(comm);
  Omega_h::test_2D_case2(comm);
}
