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

static void add_density_tag(Mesh* mesh) {
  auto density = Write<Real>(mesh->nelems());
  auto elem_coords =
      average_field(mesh, mesh->dim(), mesh->dim(), mesh->coords());
  auto f = OMEGA_H_LAMBDA(LO e) {
    if (elem_coords[e * 2 + 1] > 0.5) {
      density[e] = 1e-4;
    } else {
      density[e] = 1.0;
    }
  };
  parallel_for(mesh->nelems(), f);
  mesh->add_tag(mesh->dim(), "density", 1, Reals(density));
}

static void check_total_mass(Mesh* mesh) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  Reals masses = multiply_each(densities, sizes);
  auto owned_masses = mesh->owned_array(mesh->dim(), masses, 1);
  auto expected_mass = 1.0 * 0.5 + 1e-4 * 0.5;
  auto mass = get_sum(mesh->comm(), owned_masses);
  if (!mesh->comm()->rank()) {
    std::cout << "mass " << mass << " expected " << expected_mass << '\n';
  }
  OMEGA_H_CHECK(are_close(mass, expected_mass));
}

static Vector<2> get_total_momentum(Mesh* mesh) {
  auto vert_velocities = mesh->get_array<Real>(VERT, "velocity");
  auto elem_velocities = average_field(mesh, mesh->dim(),
      LOs(mesh->nelems(), 0, 1), mesh->dim(), vert_velocities);
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  Reals masses = multiply_each(densities, sizes);
  Reals momenta = multiply_each(elem_velocities, masses);
  auto owned_momenta = mesh->owned_array(mesh->dim(), momenta, mesh->dim());
  Vector<2> total;
  repro_sum(mesh->comm(), owned_momenta, mesh->dim(), &total[0]);
  return total;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto nx = 10;
  auto mesh = build_box(world, OMEGA_H_SIMPLEX, 1, 1, 0, nx, nx, 0);
  mesh.set_parting(OMEGA_H_GHOSTED);
  {
    auto metrics = get_implied_isos(&mesh);
    auto scalar = metric_eigenvalue_from_length(1.3);
    metrics = multiply_each_by(metrics, scalar);
    mesh.add_tag(VERT, "metric", 1, metrics);
  }
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  add_density_tag(&mesh);
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = OMEGA_H_LAMBDA(LO vert) {
    auto x = get_vector<2>(coords, vert);
    auto x2 = x[0] - 0.3;
    auto w = sign(x2) * std::sqrt(std::abs(x2));
    set_vector(velocity, vert, vector_2(1, 0) * w);
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(), Reals(velocity));
  auto momentum_before = get_total_momentum(&mesh);
  auto opts = AdaptOpts(&mesh);
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.integral_map["density"] = "mass";
  opts.xfer_opts.type_map["velocity"] = OMEGA_H_MOMENTUM_VELOCITY;
  opts.xfer_opts.velocity_density_map["velocity"] = "density";
  opts.xfer_opts.velocity_momentum_map["velocity"] = "momentum";
  opts.xfer_opts.integral_diffuse_map["mass"] =
      VarCompareOpts{VarCompareOpts::RELATIVE, 0.9, 0.0};
  opts.xfer_opts.integral_diffuse_map["momentum"] =
      VarCompareOpts{VarCompareOpts::RELATIVE, 0.02, 1e-6};
  adapt(&mesh, opts);
  check_total_mass(&mesh);
  auto momentum_after = get_total_momentum(&mesh);
  if (world->rank() == 0) {
    std::cout << "before" << ' ' << momentum_before[0] << ' '
              << momentum_before[1] << '\n';
    std::cout << "after" << ' ' << momentum_after[0] << ' ' << momentum_after[1]
              << '\n';
  }
  OMEGA_H_CHECK(are_close(momentum_before, momentum_after));
  bool ok = check_regression("gold_2d_conserve", &mesh);
  if (!ok) return 2;
  return 0;
}
