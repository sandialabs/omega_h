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

void add_dye(Mesh* mesh) {
  auto dye_w = Write<Real>(mesh->nverts());
  auto coords = mesh->coords();
  auto dye_fun = OMEGA_H_LAMBDA(LO vert) {
    auto x = get_vector<3>(coords, vert);
    auto left_cen = vector_3(.25, .5, .5);
    auto right_cen = vector_3(.75, .5, .5);
    auto left_dist = norm(x - left_cen);
    auto right_dist = norm(x - right_cen);
    auto dist = min2(left_dist, right_dist);
    if (dist < .25) {
      auto dir = sign(left_dist - right_dist);
      dye_w[vert] = 4.0 * dir * (.25 - dist);
    } else {
      dye_w[vert] = 0;
    }
  };
  parallel_for(mesh->nverts(), dye_fun);
  mesh->add_tag(VERT, "dye", 1, Reals(dye_w));
}

Reals form_pointwise(Mesh* mesh) {
  auto dim = mesh->dim();
  auto ecoords =
      average_field(mesh, dim, LOs(mesh->nelems(), 0, 1), dim, mesh->coords());
  auto pw_w = Write<Real>(mesh->nelems());
  auto pw_fun = OMEGA_H_LAMBDA(LO elem) { pw_w[elem] = ecoords[elem * dim]; };
  parallel_for(mesh->nelems(), pw_fun);
  return pw_w;
}

static void add_pointwise(Mesh* mesh) {
  auto data = form_pointwise(mesh);
  mesh->add_tag(mesh->dim(), "pointwise", 1, data);
}

static void check_total_mass(Mesh* mesh) {
  auto densities = mesh->get_array<Real>(mesh->dim(), "density");
  auto sizes = mesh->ask_sizes();
  Reals masses = multiply_each(densities, sizes);
  auto owned_masses = mesh->owned_array(mesh->dim(), masses, 1);
  OMEGA_H_CHECK(are_close(1.0, get_sum(mesh->comm(), owned_masses)));
}

static void postprocess_pointwise(Mesh* mesh) {
  auto data = mesh->get_array<Real>(mesh->dim(), "pointwise");
  auto expected = form_pointwise(mesh);
  auto diff = subtract_each(data, expected);
  mesh->add_tag(mesh->dim(), "pointwise_err", 1, diff);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  constexpr Int dim = 3;
  auto nx = 10;
  auto mesh =
      build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
  mesh.set_parting(OMEGA_H_GHOSTED);
  auto metrics = get_implied_isos(&mesh);
  mesh.add_tag(VERT, "metric", 1, metrics);
  add_dye(&mesh);
  mesh.add_tag(mesh.dim(), "density", 1, Reals(mesh.nelems(), 1.0));
  add_pointwise(&mesh);
  auto opts = AdaptOpts(&mesh);
  opts.xfer_opts.type_map["density"] = OMEGA_H_CONSERVE;
  opts.xfer_opts.integral_map["density"] = "mass";
  opts.xfer_opts.type_map["pointwise"] = OMEGA_H_POINTWISE;
  opts.xfer_opts.type_map["dye"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.integral_diffuse_map["mass"] = VarCompareOpts::none();
  opts.verbosity = EXTRA_STATS;
  auto mid = zero_vector<dim>();
  mid[0] = mid[1] = .5;
  Now t0 = now();
  for (Int i = 0; i < 8; ++i) {
    auto coords = mesh.coords();
    Write<Real> warp_w(mesh.nverts() * dim);
    auto warp_fun = OMEGA_H_LAMBDA(LO vert) {
      auto x0 = get_vector<3>(coords, vert);
      auto x1 = zero_vector<3>();
      x1[0] = x0[0];
      x1[1] = x0[1];
      auto x2 = x1 - mid;
      auto polar_a = std::atan2(x2[1], x2[0]);
      auto polar_r = norm(x2);
      Real rot_a = 0;
      if (polar_r < 0.5) {
        rot_a = (PI / 8) * (2.0 * (0.5 - polar_r));
        if (i >= 4) rot_a = -rot_a;
      }
      auto dest_a = polar_a + rot_a;
      auto dst = x0;
      dst[0] = std::cos(dest_a) * polar_r;
      dst[1] = std::sin(dest_a) * polar_r;
      dst = dst + mid;
      auto w = dst - x0;
      set_vector<3>(warp_w, vert, w);
    };
    parallel_for(mesh.nverts(), warp_fun);
    mesh.add_tag(VERT, "warp", dim, Reals(warp_w));
    while (warp_to_limit(&mesh, opts)) {
      adapt(&mesh, opts);
    }
  }
  Now t1 = now();
  mesh.set_parting(OMEGA_H_ELEM_BASED);
  if (mesh.comm()->rank() == 0) {
    std::cout << "test took " << (t1 - t0) << " seconds\n";
  }
  check_total_mass(&mesh);
  postprocess_pointwise(&mesh);
  bool ok = check_regression("gold_warp", &mesh);
  if (!ok) return 2;
  return 0;
}
