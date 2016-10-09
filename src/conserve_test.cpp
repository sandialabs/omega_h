#include <iostream>

#include "access.hpp"
#include "array.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "size.hpp"
#include "space.hpp"
#include "timer.hpp"

using namespace Omega_h;

static void postprocess_conserve(Mesh* mesh) {
  auto volume = measure_elements_real(mesh);
  auto mass = mesh->get_array<Real>(mesh->dim(), "mass");
  CHECK(are_close(1.0, sum(mesh->comm(), mass)));
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
  Vector<2> total;
  repro_sum(mesh->comm(), momenta, mesh->dim(), &total[0]);
  return total;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  CHECK(world->size() == 1);
  Mesh mesh;
  auto nx = 10;
  build_box(&mesh, lib, 1, 1, 0, nx, nx, 0);
  classify_by_angles(&mesh, PI / 4);
  mesh.reorder();
  mesh.reset_globals();
  auto size = find_implied_size(&mesh);
  size = multiply_each_by(1.3, size);
  mesh.add_tag(VERT, "size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT, size);
  mesh.add_tag(mesh.dim(), "mass", 1, OMEGA_H_CONSERVE, OMEGA_H_DO_OUTPUT,
      measure_elements_real(&mesh));
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = LAMBDA(LO vert) {
    auto x = get_vector<2>(coords, vert);
    set_vector(velocity, vert, vector_2(1, 0) * sqrt(fabs(x[0])));
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(), OMEGA_H_MOMENTUM_VELOCITY,
      OMEGA_H_DO_OUTPUT, Reals(velocity));
  auto momentum_before = get_total_momentum(&mesh);
  adapt(&mesh, AdaptOpts());
  postprocess_conserve(&mesh);
  auto momentum_after = get_total_momentum(&mesh);
  std::cout << "before" << ' ' << momentum_before[0] << ' '
            << momentum_before[1] << '\n';
  std::cout << "after" << ' ' << momentum_after[0] << ' ' << momentum_after[1]
            << '\n';
  CHECK(are_close(momentum_before, momentum_after));
  bool ok = check_regression("gold_conserve", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}
