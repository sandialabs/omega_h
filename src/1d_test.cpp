#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_internal.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_timer.hpp"

using namespace Omega_h;

static Reals logistic_function(Reals x, Real x0, Real L, Real k) {
  Write<Real> out(x.size());
  auto f = OMEGA_H_LAMBDA(LO i) {
    out[i] = L / (1 + exp(-k * (x[i] - x0)));
  };
  parallel_for(x.size(), f);
  return out;
}

static void add_solution(Mesh* mesh) {
  auto coords = mesh->coords();
  auto sol = logistic_function(coords, 0.5, 1.0, 2.0);
  mesh->add_tag(VERT, "solution", 1, sol);
}

static void add_target_metric(Mesh* mesh) {
  MetricInput input;
  input.sources.push_back({OMEGA_H_HESSIAN, true, "solution", 1.0});
  input.should_limit_element_count = true;
  input.max_element_count = 100;
  generate_target_metric_tag(mesh, input);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh(&lib);
  if (world->rank() == 0) {
    auto nx = 10;
    build_box(&mesh, 1, 0, 0, nx, 0, 0);
    classify_by_angles(&mesh, PI / 4);
    mesh.reorder();
    mesh.reset_globals();
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto metrics = get_implied_isos(&mesh);
  mesh.add_tag(VERT, "metric", 1, metrics);
  mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
  auto opts = AdaptOpts(&mesh);
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.integral_map["density"] = "mass";
  opts.xfer_opts.integral_diffuse_map["mass"] = VarCompareOpts::none();
  opts.verbosity = EXTRA_STATS;
  Now t0 = now();
  add_solution(&mesh);
  add_target_metric(&mesh);
  vtk::Writer writer(&mesh, "adapting", mesh.dim());
  writer.write();
  while (approach_metric(&mesh, opts)) {
    adapt(&mesh, opts);
    add_solution(&mesh);
    if (mesh.has_tag(VERT, "target_metric")) {
      mesh.remove_tag(VERT, "target_metric");
      add_target_metric(&mesh);
    }
    writer.write();
  }
  Now t1 = now();
  std::cout << "total time: " << (t1 - t0) << " seconds\n";
}
