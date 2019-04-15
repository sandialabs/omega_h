#include <Omega_h_amr.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <iostream>

namespace Omega_h {

static void test_2D(CommPtr comm) {
  auto mesh = build_box(comm, OMEGA_H_HYPERCUBE,
      2, 1., 0., 2, 1, 0); 
  vtk::FullWriter writer("out_amr_mpi_2D", &mesh);
  if (!comm->rank()) std::cout << "before first adapt" << std::endl;
  writer.write();
  { // first adapt
    Write<Byte> marks(mesh.nelems(), 0); 
    if (comm->rank() == 1) marks.set(0, 1); 
    const auto xfer_opts = TransferOpts();
    amr::refine(&mesh, marks, xfer_opts);
  }
  if (!comm->rank()) std::cout << "before second adapt" << std::endl;
  writer.write();
  { // second adapt
    Write<Byte> marks(mesh.nelems(), 0); 
    if (comm->rank() == 1) marks.set(2, 1);
    const auto xfer_opts = TransferOpts();
    amr::refine(&mesh, marks, xfer_opts);
  }
}

}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  OMEGA_H_CHECK(comm->size() == 2); 
  Omega_h::test_2D(comm);
}
