#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_print.hpp>
#include <Omega_h_file.hpp>

#include <Omega_h_print.hpp>

struct MyWriter {
  MyWriter(Omega_h::Mesh* m)
      : face("outamr_face", m, 2),
        edge("outamr_edge", m, 1) {}
  Omega_h::vtk::Writer face;
  Omega_h::vtk::Writer edge;
  void write() {
    face.write();
    edge.write();
  }
};

static void write_parent(
    Omega_h::LOs parent_idx,
    Omega_h::Read<Omega_h::I8> codes) {
  OMEGA_H_CHECK(parent_idx.size() == codes.size());
  std::cout << "i        parent_idx      which_child     parent_dim" << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  for (Omega_h::LO i = 0; i < parent_idx.size(); ++i) {
    std::cout << i << "\t\t" << parent_idx[i] << "\t\t"
      << Omega_h::code_which_child(codes[i]) << "\t\t"
      << Omega_h::code_parent_dim(codes[i]) << std::endl;
  }
}

static void write_parents(Omega_h::Mesh* mesh) {
  for (Omega_h::Int i = 0; i <= mesh->dim(); ++i) {
    std::cout << "child dim:" << i << std::endl;
    std::cout << "-------------" << std::endl;
    auto parents = mesh->ask_parents(i);
    auto parent_idx = parents.parent_idx;
    auto codes = parents.codes;
    write_parent(parent_idx, codes);
    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto mesh = Omega_h::build_box(
      lib.world(), OMEGA_H_HYPERCUBE, 2.0, 1.0, 0.0, 2, 1, 0);
  MyWriter writer(&mesh);
  writer.write();
  {
    std::cout << "adapt 1" << std::endl;
    auto xfer_opts = Omega_h::TransferOpts();
    Omega_h::amr_refine(&mesh, Omega_h::Bytes({1, 0}), xfer_opts);
    writer.write();
    write_parents(&mesh);
  }
  {
    std::cout << "adapt 1" << std::endl;
    auto xfer_opts = Omega_h::TransferOpts();
    Omega_h::amr_refine(&mesh, Omega_h::Bytes({0,1,0,0,0,0}), xfer_opts);
    writer.write();
    write_parents(&mesh);
  }
  {
    std::cout << "adapt 1" << std::endl;
    auto xfer_opts = Omega_h::TransferOpts();
    Omega_h::amr_refine(&mesh, Omega_h::Bytes({0,0,0,0,0,1,0,0,1,0}), xfer_opts);
    writer.write();
    write_parents(&mesh);
  }
}
