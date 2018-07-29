#include <iomanip>
#include <iostream>

#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_timer.hpp"

using namespace Omega_h;

constexpr Int nobjs = 2;
constexpr Int obj_ids[nobjs] = {34, 72};

static void check_total_mass(Mesh* mesh) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  Reals masses = multiply_each(densities, sizes);
  auto owned_masses = mesh->owned_array(mesh->dim(), masses, 1);
  auto total_mass = get_sum(mesh->comm(), owned_masses);
  if (mesh->comm()->rank() == 0) {
    std::cerr << "total mass " << total_mass << '\n';
  }
  OMEGA_H_CHECK(are_close(1.0, total_mass, 1e-4, 1e-4));
}

static Real get_object_size(Mesh* mesh, Int obj) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  auto masses = multiply_each(densities, sizes);
  auto class_ids = mesh->get_array<ClassId>(mesh->dim(), "class_id");
  auto elem_in_obj = each_eq_to(class_ids, obj_ids[obj]);
  auto obj_elems = collect_marked(elem_in_obj);
  auto obj_sizes = unmap(obj_elems, sizes, 1);
  return repro_sum(mesh->comm(), obj_sizes);
}

static Real get_object_mass(Mesh* mesh, Int obj) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  Reals masses = multiply_each(densities, sizes);
  auto class_ids = mesh->get_array<ClassId>(mesh->dim(), "class_id");
  auto elem_in_obj = each_eq_to(class_ids, obj_ids[obj]);
  auto obj_elems = collect_marked(elem_in_obj);
  auto obj_masses = unmap(obj_elems, masses, 1);
  return repro_sum(mesh->comm(), obj_masses);
}

static Vector<3> get_total_momentum(Mesh* mesh) {
  auto vert_velocities = mesh->get_array<Real>(VERT, "velocity");
  auto elem_velocities =
      average_field(mesh, mesh->dim(), mesh->dim(), vert_velocities);
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  Reals masses = multiply_each(densities, sizes);
  Reals momenta = multiply_each(elem_velocities, masses);
  Vector<3> total;
  repro_sum(mesh->comm(), momenta, mesh->dim(), &total[0]);
  return total;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 2);
  auto world = lib.world();
  auto mesh = gmsh::read(argv[1], world);
  mesh.set_parting(OMEGA_H_GHOSTED);
  {
    auto metrics = get_implied_isos(&mesh);
    auto scalar = metric_eigenvalue_from_length(1.2);
    metrics = multiply_each_by(metrics, scalar);
    mesh.add_tag(VERT, "metric", 1, metrics);
  }
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
  auto velocity = Write<Real>(mesh.nverts() * mesh.dim());
  auto coords = mesh.coords();
  auto f = OMEGA_H_LAMBDA(LO vert) {
    auto x = get_vector<3>(coords, vert);
    set_vector(velocity, vert, vector_3(0, 0, 1) * std::sqrt(std::abs(x[2])));
  };
  parallel_for(mesh.nverts(), f);
  mesh.add_tag(VERT, "velocity", mesh.dim(), Reals(velocity));
  fix_momentum_velocity_verts(&mesh, {{2, 10}}, 2);
  auto momentum_before = get_total_momentum(&mesh);
  Real masses_before[nobjs];
  Real sizes_before[nobjs];
  for (Int obj = 0; obj < nobjs; ++obj) {
    masses_before[obj] = get_object_mass(&mesh, obj);
    sizes_before[obj] = get_object_size(&mesh, obj);
  }
  auto opts = AdaptOpts(&mesh);
  opts.verbosity = EXTRA_STATS;
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.type_map["velocity"] = OMEGA_H_MOMENTUM_VELOCITY;
  opts.xfer_opts.integral_map["density"] = "mass";
  opts.xfer_opts.velocity_density_map["velocity"] = "density";
  opts.xfer_opts.velocity_momentum_map["velocity"] = "momentum";
  opts.xfer_opts.integral_diffuse_map["mass"] =
      VarCompareOpts{VarCompareOpts::RELATIVE, 0.2, 0.0};
  opts.xfer_opts.integral_diffuse_map["momentum"] =
      VarCompareOpts{VarCompareOpts::RELATIVE, 0.02, 1e-6};
  opts.xfer_opts.should_conserve_size = true;
  adapt(&mesh, opts);
  check_total_mass(&mesh);
  for (Int obj = 0; obj < nobjs; ++obj) {
    auto mass_after = get_object_mass(&mesh, obj);
    auto size_after = get_object_size(&mesh, obj);
    if (world->rank() == 0) {
      std::cout << "model region " << obj_ids[obj] << " mass before "
                << masses_before[obj] << ", after " << mass_after << '\n';
      std::cout << "model region " << obj_ids[obj] << " size before "
                << sizes_before[obj] << ", after " << size_after << '\n';
    }
    OMEGA_H_CHECK(are_close(mass_after, masses_before[obj], 1e-4, 1e-4));
  }
  auto momentum_after = get_total_momentum(&mesh);
  if (world->rank() == 0) {
    std::cout << "momentum before" << ' ' << momentum_before[0] << ' '
              << momentum_before[1] << ' ' << momentum_before[2] << '\n';
    std::cout << "momentum after" << ' ' << momentum_after[0] << ' '
              << momentum_after[1] << ' ' << momentum_after[2] << '\n';
  }
  OMEGA_H_CHECK(are_close(momentum_before, momentum_after));
  bool ok = check_regression("gold_3d_conserve", &mesh);
  if (!ok) return 2;
  return 0;
}
