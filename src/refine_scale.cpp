#include <Omega_h_build.hpp>
#include <Omega_h_adapt.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_stack.hpp>

#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::ScopedTimer scoped_timer("main");
  auto world = lib.world();
  if (world->rank() == 0) std::cout << "building a box of " << (6 * 7 * 7 * 210) << " elements!\n";
  auto mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1e-6, 1e-6, 1e-3, 7, 7, 210);
  if (world->rank() == 0) std::cout << "done building the box, re-balancing!\n";
  mesh.balance();
  if (world->rank() == 0) std::cout << "done rebalancing, rank 0 has " << mesh.nelems() << " elements!\n";
  if (world->rank() == 0) std::cout << "ghosting!\n";
  mesh.set_parting(OMEGA_H_GHOSTED);
  if (world->rank() == 0) std::cout << "done ghosting, rank 0 has " << mesh.nelems() << " elements!\n";
  if (world->rank() == 0) std::cout << "computing implied metric!\n";
  Omega_h::add_implied_metric_tag(&mesh);
  auto ncurrent_elems = mesh.nglobal_ents(mesh.dim());
  auto ntarget_elems = ncurrent_elems * 8;
  if (world->rank() == 0) std::cout << "done with implied metric, scaling it to " << ntarget_elems << " elements!\n"; 
  auto metrics = mesh.get_array<double>(0, "metric");
  auto scalar = Omega_h::get_metric_scalar_for_nelems(mesh.dim(), ncurrent_elems, ntarget_elems);
  if (world->rank() == 0) std::cout << "scalar is " << scalar << "!\n";
  metrics = multiply_each_by(metrics, scalar);
  mesh.add_tag(Omega_h::VERT, "target_metric", Omega_h::symm_ncomps(mesh.dim()), metrics);
  auto opts = Omega_h::AdaptOpts(&mesh);
  opts.verbosity = Omega_h::EXTRA_STATS;
  if (world->rank() == 0) std::cout << "done scaling, adapting!\n";
  while (Omega_h::approach_metric(&mesh, opts)) {
      Omega_h::adapt(&mesh, opts);
  }
}
