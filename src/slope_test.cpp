#include "Omega_h.hpp"
#include "Omega_h_math.hpp"
#include "Omega_h_array_ops.hpp"
#include "loop.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  char const* path_in = nullptr;
  char const* path_out = nullptr;
  char const* vtk_path = nullptr;
  bool usage = false;
  for (int i = 1; i < argc; ++i) {
    if (!strcmp("-g", argv[i])) {
      if (i == argc - 1) {
        std::cout << " -g takes an argument\n";
        usage = true;
      }
      vtk_path = argv[++i];
    } else if (!path_in)
      path_in = argv[i];
    else if (!path_out)
      path_out = argv[i];
  }
  if (!path_in || !path_out) usage = true;
  if (usage) {
    std::cout << "usage: " << argv[0] << " [options] input.osh output.osh\n";
    std::cout << "  -g <out_vtk>\n";
    return -1;
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(path_in, lib.world(), &mesh);
  auto target_nelems =
      mesh.comm()->allreduce<Omega_h::Real>(mesh.nelems(), OMEGA_H_SUM);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto dim = mesh.dim();
  auto coords = mesh.coords();
  auto bb = Omega_h::get_bounding_box<3>(&mesh);
  auto analytic_size_w = Omega_h::Write<Omega_h::Real>(mesh.nverts());
  auto f = LAMBDA(Omega_h::LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto minz = bb.min[dim - 1];
    auto maxz = bb.max[dim - 1];
    auto z_norm = (z * 2.0 / (maxz - minz)) - 1.0;
    auto bot_h = -2.0 * z_norm + 1.0;
    auto top_h = 2.0 * z_norm + 1.0;
    auto h = Omega_h::max2(bot_h, top_h);
    analytic_size_w[v] = h;
  };
  Omega_h::parallel_for(mesh.nverts(), f);
  auto analytic_size = Omega_h::Reals(analytic_size_w);
  Omega_h::vtk::Writer writer;
  if (vtk_path) writer = Omega_h::vtk::Writer(&mesh, vtk_path, dim);
  mesh.add_tag(
      Omega_h::VERT, "size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT, analytic_size);
  if (vtk_path) writer.write();
  auto scalar =
      Omega_h::size_scalar_for_nelems(&mesh, analytic_size, target_nelems);
  auto scaled_size = Omega_h::multiply_each_by(scalar, analytic_size);
  mesh.set_tag(Omega_h::VERT, "size", scaled_size);
  auto imb = mesh.imbalance();
  if (!mesh.comm()->rank()) std::cout << "imbalance on input " << imb << '\n';
  if (vtk_path) writer.write();
  mesh.balance(true);
  imb = mesh.imbalance();
  if (!mesh.comm()->rank())
    std::cout << "imbalance after predictive " << imb << '\n';
  if (vtk_path) writer.write();
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;
  Omega_h::adapt(&mesh, opts);
  imb = mesh.imbalance();
  if (!mesh.comm()->rank())
    std::cout << "imbalance after adapt " << imb << '\n';
  if (vtk_path) writer.write();
  mesh.balance();
  imb = mesh.imbalance();
  if (!mesh.comm()->rank())
    std::cout << "imbalance after post-balance " << imb << '\n';
  if (vtk_path) writer.write();
  mesh.remove_tag(Omega_h::VERT, "size");
  Omega_h::binary::write(path_out, &mesh);
}
