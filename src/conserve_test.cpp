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

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  Mesh mesh;
  if (world->rank() == 0) {
    auto nx = 10;
    build_box(&mesh, lib, 1, 1, 0, nx, nx, 0);
    classify_by_angles(&mesh, PI / 4);
    mesh.reorder();
    mesh.reset_globals();
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto size = find_identity_size(&mesh);
  size = multiply_each_by(2.0, size);
  mesh.add_tag(VERT, "size", 1, OMEGA_H_LINEAR_INTERP, OMEGA_H_DO_OUTPUT, size);
  mesh.add_tag(mesh.dim(), "mass", 1, OMEGA_H_CONSERVE, OMEGA_H_DO_OUTPUT,
      measure_elements_real(&mesh));
  mesh.add_tag(VERT, "velocity", mesh.dim(), OMEGA_H_MOMENTUM_VELOCITY,
    OMEGA_H_DO_OUTPUT, Reals(mesh.nverts() * mesh.dim(), 1));
  adapt(&mesh, 0.30, 0.30, 2.0 / 3.0, 4.0 / 3.0, 4, 3);
  postprocess_conserve(&mesh);
  bool ok = check_regression("gold_conserve", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}

