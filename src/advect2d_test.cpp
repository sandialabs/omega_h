#include "Omega_h.hpp"
#include "Omega_h_compare.hpp"

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 2);
  auto inpath = argv[1];
  auto mesh = Omega_h::Mesh(&lib);
  auto maximum_size = Omega_h::Real(0.9);
  auto target_error = Omega_h::Real(0.011);
  auto gradation_rate = Omega_h::Real(1.0);
  auto max_metric_length = Omega_h::Real(2.8);
  Omega_h::binary::read(inpath, lib.world(), &mesh);
  mesh.balance();
  mesh.reorder();
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto sol = mesh.get_array<Omega_h::Real>(Omega_h::VERT, "Solution");
  auto hessians = Omega_h::recover_hessians(&mesh, sol);
  auto metrics =
      Omega_h::metric_from_hessians(mesh.dim(), hessians, target_error);
  metrics = Omega_h::clamp_metrics(mesh.nverts(), metrics, 0.0, maximum_size);
  metrics = Omega_h::limit_metric_gradation(&mesh, metrics, gradation_rate);
  mesh.add_tag(
      Omega_h::VERT, "target_metric", Omega_h::symm_ncomps(mesh.dim()), metrics);
  auto implied_metrics = Omega_h::find_implied_metric(&mesh);
  mesh.add_tag(
      Omega_h::VERT, "metric", Omega_h::symm_ncomps(mesh.dim()), implied_metrics);
  Omega_h::AdaptOpts opts(&mesh);
  opts.xfer_opts.type_map["Solution"] = OMEGA_H_LINEAR_INTERP;
  opts.max_length_allowed = max_metric_length;
  while (Omega_h::approach_metric(&mesh, opts)) {
    Omega_h::adapt(&mesh, opts);
  }
  bool ok = check_regression("gold_advect2d", &mesh);
  if (!ok) return 2;
  return 0;
}
