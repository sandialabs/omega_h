#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_print.hpp>
#include <Omega_h_file.hpp>

int main(int argc, char** argv) {

  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(
      lib.world(), OMEGA_H_HYPERCUBE, 2.0, 1.0, 0.0, 2, 1, 0);
  auto xfer_opts = Omega_h::TransferOpts();
  Omega_h::amr_refine(&mesh, Omega_h::Bytes({1, 0}), xfer_opts);

  auto face_writer = Omega_h::vtk::Writer("face_outamr", &mesh, 2);
  auto edge_writer = Omega_h::vtk::Writer("edge_outamr", &mesh, 1);
  auto vtx_writer = Omega_h::vtk::Writer("vtx_outamr", &mesh, 0);

  face_writer.write();
  edge_writer.write();
  vtx_writer.write();

  Omega_h::amr_refine(&mesh, Omega_h::Bytes({0,1,0,0,0,0}), xfer_opts);

  face_writer.write();
  edge_writer.write();
  vtx_writer.write();

  Omega_h::amr_refine(&mesh, Omega_h::Bytes({0,0,0,0,0,1,0,0,1,0}), xfer_opts);

  face_writer.write();
  edge_writer.write();
  vtx_writer.write();

}
