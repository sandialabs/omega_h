#include <Omega_h_adapt.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_timer.hpp>

#include <iostream>

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto mesh = Mesh(&lib);
  binary::read(argv[1], lib.world(), &mesh);
  auto all_h = mesh.get_array<Real>(VERT, "proteus_size_scale");
  auto all_R = mesh.get_array<Real>(VERT, "proteus_size_frame");
  auto all_M_without_frames =
      Write<Real>(mesh.nverts() * symm_ncomps(mesh.dim()));
  auto all_M_with_frames = Write<Real>(mesh.nverts() * symm_ncomps(mesh.dim()));
  auto all_M_without_scales =
      Write<Real>(mesh.nverts() * symm_ncomps(mesh.dim()));
  constexpr auto dim = 3;
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto h = get_vector<dim>(all_h, v);
    Matrix<dim, dim> R;
    for (Int j = 0; j < dim; ++j) {
      for (Int k = 0; k < dim; ++k) {
        R[j][k] = all_R[v * dim * dim + k * dim + j];
      }
    }
    auto M_without_frame = diagonal(metric_eigenvalues_from_lengths(h));
    auto M_with_frame = compose_metric(R, h);
    auto M_without_scale = R * transpose(R);
    set_symm(all_M_without_frames, v, M_without_frame);
    set_symm(all_M_with_frames, v, M_with_frame);
    set_symm(all_M_without_scales, v, M_without_scale);
  };
  parallel_for(mesh.nverts(), f, "convert_metric");
  mesh.add_tag(
      VERT, "M_without_frame", symm_ncomps(dim), Reals(all_M_without_frames));
  mesh.add_tag(
      VERT, "M_with_frame", symm_ncomps(dim), Reals(all_M_with_frames));
  mesh.add_tag(
      VERT, "M_without_scales", symm_ncomps(dim), Reals(all_M_without_scales));
  vtk::write_vtu("converted.vtu", &mesh);
  mesh.remove_tag(VERT, "M_without_frame");
  mesh.remove_tag(VERT, "M_with_frame");
  mesh.remove_tag(VERT, "M_without_scales");
#if 1
  auto target_metric = Reals(all_M_without_frames);
#else
  auto target_metric = Reals(all_M_with_frames);
#endif
  mesh.add_tag(VERT, "original_metric", symm_ncomps(dim), target_metric);
  auto t0 = now();
  target_metric = limit_metric_gradation(&mesh, target_metric, 1.0);
  mesh.add_tag(VERT, "target_metric", symm_ncomps(dim), target_metric);
  add_implied_metric_tag(&mesh);
  auto t1 = now();
  auto opts = AdaptOpts(&mesh);
  vtk::Writer writer("adapting", &mesh);
  mesh.ask_qualities();
  writer.write();
  auto t2 = now();
  fix(&mesh, opts, OMEGA_H_ANISOTROPIC, true);
  auto t3 = now();
  Real adapt_time = 0;
  while (true) {
    auto t4 = now();
    if (!approach_metric(&mesh, opts)) break;
    adapt(&mesh, opts);
    auto t5 = now();
    adapt_time += (t5 - t4);
    writer.write();
  }
  std::cout << "size field preprocessing time: " << (t1 - t0) << " seconds\n";
  std::cout << "initial quality repair time: " << (t3 - t2) << " seconds\n";
  std::cout << "approach and adapt time: " << adapt_time << " seconds\n";
}
