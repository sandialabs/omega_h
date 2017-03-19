#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_timer.hpp"
#include "access.hpp"
#include "internal.hpp"
#include "size.hpp"
#include "space.hpp"

using namespace Omega_h;

constexpr Int nobjs = 2;
constexpr Int obj_ids[nobjs] = {34, 72};

static void postprocess_conserve(Mesh* mesh) {
  auto volume = measure_elements_real(mesh);
  auto mass = mesh->get_array<Real>(mesh->dim(), "mass");
  auto owned_mass = mesh->owned_array(mesh->dim(), mass, 1);
  CHECK(are_close(1.0, get_sum(mesh->comm(), owned_mass)));
  auto density = divide_each(mass, volume);
  mesh->add_tag(mesh->dim(), "density", 1,
      density);
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
  CHECK(argc == 2);
  auto world = lib.world();
  Mesh mesh(&lib);
  if (world->rank() == 0) {
    gmsh::read(argv[1], &mesh);
    mesh.reorder();
    mesh.reset_globals();
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OMEGA_H_GHOSTED);
  {
    auto size = find_implied_size(&mesh);
    size = multiply_each_by(1.2, size);
    mesh.add_tag(VERT, "size", 1, size);
  }
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  mesh.add_tag(mesh.dim(), "mass", 1,
      measure_elements_real(&mesh));
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = LAMBDA(LO vert) {
    auto x = get_vector<3>(coords, vert);
    set_vector(velocity, vert, vector_3(0, 0, 1) * sqrt(fabs(x[2])));
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(),
      Reals(velocity));
  fix_momentum_velocity_verts(&mesh, 2, 10, 2);
  auto momentum_before = get_total_momentum(&mesh);
  Real masses_before[nobjs];
  for (Int obj = 0; obj < nobjs; ++obj) {
    masses_before[obj] = get_total_mass(&mesh, obj);
  }
  adapt(&mesh, AdaptOpts(&mesh));
  postprocess_conserve(&mesh);
  for (Int obj = 0; obj < nobjs; ++obj) {
    auto mass_after = get_total_mass(&mesh, obj);
    if (world->rank() == 0) {
      std::cout << "model region " << obj_ids[obj] << " mass before "
                << masses_before[obj] << ", after " << mass_after << '\n';
    }
    CHECK(are_close(mass_after, masses_before[obj]));
  }
  auto momentum_after = get_total_momentum(&mesh);
  if (world->rank() == 0) {
    std::cout << "momentum before" << ' ' << momentum_before[0] << ' '
              << momentum_before[1] << ' ' << momentum_before[2] << '\n';
    std::cout << "momentum after" << ' ' << momentum_after[0] << ' '
              << momentum_after[1] << ' ' << momentum_after[2] << '\n';
  }
  CHECK(are_close(momentum_before, momentum_after));
  bool ok = check_regression("gold_3d_conserve", &mesh);
  if (!ok) return 2;
  return 0;
}
