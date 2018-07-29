#include <iostream>

#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_timer.hpp"

using namespace Omega_h;

static Reals logistic_function(Reals x, Real x0, Real L, Real k) {
  Write<Real> out(x.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    out[i] = L / (1 + std::exp(-k * (x[i] - x0)));
  };
  parallel_for(x.size(), f);
  return out;
}

static void add_solution(Mesh* mesh) {
  auto coords = mesh->coords();
  auto sol = logistic_function(coords, 0.5, 1.0, 20.0);
  mesh->add_tag(VERT, "solution", 1, sol);
}

static void add_metric(Mesh* mesh) {
  MetricInput input;
  input.sources.push_back(MetricSource{OMEGA_H_VARIATION, 1.0, "solution"});
  input.should_limit_lengths = true;
  input.max_length = 1.0;
  input.should_limit_gradation = true;
  input.max_gradation_rate = 1.0;
  input.should_limit_element_count = true;
  input.max_element_count = 100;
  input.min_element_count = 50;
  generate_metric_tag(mesh, input);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto nx = 10;
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, 1., 0., 0., nx, 0, 0);
  mesh.set_parting(OMEGA_H_GHOSTED);
  mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
  auto opts = AdaptOpts(&mesh);
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.integral_map["density"] = "mass";
  opts.xfer_opts.integral_diffuse_map["mass"] = VarCompareOpts::none();
  opts.verbosity = EXTRA_STATS;
  opts.nquality_histogram_bins = 1;
  Now t0 = now();
  add_solution(&mesh);
  add_metric(&mesh);
  while (1) {
    adapt(&mesh, opts);
    mesh.set_parting(OMEGA_H_GHOSTED);
    add_solution(&mesh);
    mesh.remove_tag(VERT, "metric");
    add_metric(&mesh);
    if (mesh.max_length() < 2.0) break;
  }
  Now t1 = now();
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  if (world->rank() == 0) {
    std::cout << "total time: " << (t1 - t0) << " seconds\n";
  }
  bool ok = check_regression("gold_1d", &mesh);
  if (!ok) return 2;
  return 0;
}
