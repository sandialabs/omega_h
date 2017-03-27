#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_timer.hpp"
#include "internal.hpp"
#include "Omega_h_shape.hpp"

using namespace Omega_h;

static void check_total_mass(Mesh* mesh) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  auto masses = multiply_each(densities, sizes);
  auto owned_masses = mesh->owned_array(mesh->dim(), masses, 1);
  CHECK(are_close(1.0, get_sum(mesh->comm(), owned_masses)));
}

static Vector<2> get_total_momentum(Mesh* mesh) {
  auto vert_velocities = mesh->get_array<Real>(VERT, "velocity");
  auto elem_velocities = average_field(mesh, mesh->dim(),
      LOs(mesh->nelems(), 0, 1), mesh->dim(), vert_velocities);
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  auto masses = multiply_each(densities, sizes);
  auto momenta = multiply_each(elem_velocities, masses);
  auto owned_momenta = mesh->owned_array(mesh->dim(), momenta, mesh->dim());
  Vector<2> total;
  repro_sum(mesh->comm(), owned_momenta, mesh->dim(), &total[0]);
  return total;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh(&lib);
  if (world->rank() == 0) {
    auto nx = 10;
    build_box(&mesh, 1, 1, 0, nx, nx, 0);
    classify_by_angles(&mesh, PI / 4);
    mesh.reorder();
    mesh.reset_globals();
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  {
    auto metrics = find_implied_isos(&mesh);
    auto scalar = metric_eigenvalue_from_length(1.3);
    metrics = multiply_each_by(scalar, metrics);
    mesh.add_tag(VERT, "metric", 1, metrics);
  }
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = LAMBDA(LO vert) {
    auto x = get_vector<2>(coords, vert);
    auto x2 = x[0] - 0.3;
    auto w = sign(x2) * sqrt(fabs(x2));
    set_vector(velocity, vert, vector_2(1, 0) * w);
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(), Reals(velocity));
  auto momentum_before = get_total_momentum(&mesh);
  auto opts = AdaptOpts(&mesh);
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.type_map["velocity"] = OMEGA_H_MOMENTUM_VELOCITY;
  opts.xfer_opts.velocity_density_map["velocity"] = "density";
  opts.xfer_opts.velocity_momentum_map["velocity"] = "momentum";
  opts.xfer_opts.integral_diffuse_map["density"] = VarCompareOpts::none();
  opts.xfer_opts.integral_diffuse_map["momentum"] = VarCompareOpts::none();
  adapt(&mesh, opts);
  check_total_mass(&mesh);
  auto momentum_after = get_total_momentum(&mesh);
  if (world->rank() == 0) {
    std::cout << "before" << ' ' << momentum_before[0] << ' '
              << momentum_before[1] << '\n';
    std::cout << "after" << ' ' << momentum_after[0] << ' ' << momentum_after[1]
              << '\n';
  }
  CHECK(are_close(momentum_before, momentum_after));
  bool ok = check_regression("gold_2d_conserve", &mesh);
  if (!ok) return 2;
  return 0;
}
