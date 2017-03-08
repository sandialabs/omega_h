#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "access.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "size.hpp"
#include "space.hpp"
#include "timer.hpp"

using namespace Omega_h;

static void postprocess_conserve(Mesh* mesh) {
  auto volume = measure_elements_real(mesh);
  auto mass = mesh->get_array<Real>(mesh->dim(), "mass");
  auto owned_mass = mesh->owned_array(mesh->dim(), mass, 1);
  CHECK(are_close(1.0, get_sum(mesh->comm(), owned_mass)));
  auto density = divide_each(mass, volume);
  mesh->add_tag(mesh->dim(), "density", 1, OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DO_OUTPUT, density);
}

static Vector<2> get_total_momentum(Mesh* mesh) {
  auto vert_velocities = mesh->get_array<Real>(VERT, "velocity");
  auto elem_velocities = average_field(mesh, mesh->dim(),
      LOs(mesh->nelems(), 0, 1), mesh->dim(), vert_velocities);
  auto masses = mesh->get_array<Real>(mesh->dim(), "mass");
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
    auto size = find_implied_size(&mesh);
    size = multiply_each_by(1.3, size);
    mesh.add_tag(VERT, "size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT, size);
  }
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  mesh.add_tag(mesh.dim(), "mass", 1, OMEGA_H_CONSERVE, OMEGA_H_DO_OUTPUT,
      measure_elements_real(&mesh));
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = LAMBDA(LO vert) {
    auto x = get_vector<2>(coords, vert);
    auto x2 = x[0] - 0.3;
    auto w = sign(x2) * sqrt(fabs(x2));
    set_vector(velocity, vert, vector_2(1, 0) * w);
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(), OMEGA_H_MOMENTUM_VELOCITY,
      OMEGA_H_DO_OUTPUT, Reals(velocity));
  auto momentum_before = get_total_momentum(&mesh);
  adapt(&mesh, AdaptOpts(&mesh));
  postprocess_conserve(&mesh);
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
