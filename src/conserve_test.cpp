#include <iostream>

#include "access.hpp"
#include "array.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "map.hpp"
#include "size.hpp"
#include "space.hpp"
#include "timer.hpp"

using namespace Omega_h;

constexpr Int nobjs = 2;
constexpr Int obj_ids[nobjs] = {34, 72};

static void postprocess_conserve(Mesh* mesh) {
  auto volume = measure_elements_real(mesh);
  auto mass = mesh->get_array<Real>(mesh->dim(), "mass");
  CHECK(are_close(1.0, sum(mesh->comm(), mass)));
  auto density = divide_each(mass, volume);
  mesh->add_tag(mesh->dim(), "density", 1, OMEGA_H_DONT_TRANSFER,
      OMEGA_H_DO_OUTPUT, density);
}

static Real get_total_mass(Mesh* mesh, Int obj) {
  auto masses = mesh->get_array<Real>(mesh->dim(), "mass");
  auto class_ids = mesh->get_array<I32>(mesh->dim(), "class_id");
  auto elem_in_obj = each_eq_to(class_ids, obj_ids[obj]);
  auto obj_elems = collect_marked(elem_in_obj);
  auto obj_masses = unmap(obj_elems, masses, 1);
  return repro_sum(mesh->comm(), obj_masses);
}

static Vector<3> get_total_momentum(Mesh* mesh) {
  auto vert_velocities = mesh->get_array<Real>(VERT, "velocity");
  auto elem_velocities = average_field(mesh, mesh->dim(),
      LOs(mesh->nelems(), 0, 1), mesh->dim(), vert_velocities);
  auto masses = mesh->get_array<Real>(mesh->dim(), "mass");
  auto momenta = multiply_each(elem_velocities, masses);
  Vector<3> total;
  repro_sum(mesh->comm(), momenta, mesh->dim(), &total[0]);
  return total;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  CHECK(world->size() == 1);
  Mesh mesh(&lib);
  gmsh::read("ball_in_cube.msh", &mesh);
  auto size = find_implied_size(&mesh);
  size = multiply_each_by(1.2, size);
  mesh.add_tag(VERT, "size", 1, OMEGA_H_SIZE, OMEGA_H_DO_OUTPUT, size);
  mesh.add_tag(mesh.dim(), "mass", 1, OMEGA_H_CONSERVE, OMEGA_H_DO_OUTPUT,
      measure_elements_real(&mesh));
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = LAMBDA(LO vert) {
    auto x = get_vector<3>(coords, vert);
    set_vector(velocity, vert, vector_3(1, 0, 0) * sqrt(fabs(x[0])));
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(), OMEGA_H_MOMENTUM_VELOCITY,
      OMEGA_H_DO_OUTPUT, Reals(velocity));
//fix_momentum_velocity_verts(
//    &mesh, std::vector<Int>({2}), std::vector<I32>({33}));
  auto momentum_before = get_total_momentum(&mesh);
  Real masses_before[nobjs];
  for (Int obj = 0; obj < nobjs; ++obj) {
    masses_before[obj] = get_total_mass(&mesh, obj);
  }
  adapt(&mesh, AdaptOpts(&mesh));
  postprocess_conserve(&mesh);
  for (Int obj = 0; obj < nobjs; ++obj) {
    auto mass_after = get_total_mass(&mesh, obj);
    std::cout << "model region " << obj_ids[obj] << " mass before "
              << masses_before[obj] << ", after " << mass_after << '\n';
    CHECK(are_close(mass_after, masses_before[obj]));
  }
  auto momentum_after = get_total_momentum(&mesh);
  std::cout << "momentum before" << ' ' << momentum_before[0] << ' '
            << momentum_before[1] << ' ' << momentum_before[2] << '\n';
  std::cout << "momentum after" << ' ' << momentum_after[0] << ' '
            << momentum_after[1] << ' ' << momentum_after[2] << '\n';
  CHECK(are_close(momentum_before, momentum_after));
  bool ok = check_regression("gold_conserve", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}
