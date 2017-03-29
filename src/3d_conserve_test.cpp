#include <iostream>
#include <iomanip>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_timer.hpp"
#include "internal.hpp"
#include "Omega_h_shape.hpp"

using namespace Omega_h;

constexpr Int nobjs = 2;
constexpr Int obj_ids[nobjs] = {34, 72};

static void check_total_mass(Mesh* mesh) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  auto masses = multiply_each(densities, sizes);
  auto owned_masses = mesh->owned_array(mesh->dim(), masses, 1);
  OMEGA_H_CHECK(are_close(1.0, get_sum(mesh->comm(), owned_masses)));
}

static Real get_object_mass(Mesh* mesh, Int obj) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  auto masses = multiply_each(densities, sizes);
  auto class_ids = mesh->get_array<I32>(mesh->dim(), "class_id");
  auto elem_in_obj = each_eq_to(class_ids, obj_ids[obj]);
  auto obj_elems = collect_marked(elem_in_obj);
  auto obj_masses = unmap(obj_elems, masses, 1);
  return repro_sum(mesh->comm(), obj_masses);
}

static Vector<3> get_total_momentum(Mesh* mesh) {
  auto vert_velocities = mesh->get_array<Real>(VERT, "velocity");
  auto elem_velocities = average_field(
      mesh, mesh->dim(), mesh->dim(), vert_velocities);
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  auto masses = multiply_each(densities, sizes);
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
    auto metrics = find_implied_isos(&mesh);
    auto scalar = metric_eigenvalue_from_length(1.2);
    metrics = multiply_each_by(scalar, metrics);
    mesh.add_tag(VERT, "metric", 1, metrics);
  }
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = LAMBDA(LO vert) {
    auto x = get_vector<3>(coords, vert);
    set_vector(velocity, vert, vector_3(0, 0, 1) * sqrt(fabs(x[2])));
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(), Reals(velocity));
  fix_momentum_velocity_verts(&mesh, 2, 10, 2);
  auto momentum_before = get_total_momentum(&mesh);
  Real masses_before[nobjs];
  for (Int obj = 0; obj < nobjs; ++obj) {
    masses_before[obj] = get_object_mass(&mesh, obj);
  }
  auto opts = AdaptOpts(&mesh);
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.type_map["velocity"] = OMEGA_H_MOMENTUM_VELOCITY;
  opts.xfer_opts.integral_map["density"] = "mass";
  opts.xfer_opts.velocity_density_map["velocity"] = "density";
  opts.xfer_opts.velocity_momentum_map["velocity"] = "momentum";
  opts.xfer_opts.integral_diffuse_map["mass"] = VarCompareOpts::none();
  opts.xfer_opts.integral_diffuse_map["momentum"] = VarCompareOpts::none();
  adapt(&mesh, opts);
  check_total_mass(&mesh);
  for (Int obj = 0; obj < nobjs; ++obj) {
    auto mass_after = get_object_mass(&mesh, obj);
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
