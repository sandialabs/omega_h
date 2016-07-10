#include <iostream>

#include "access.hpp"
#include "array.hpp"
#include "classify.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "size.hpp"
#include "space.hpp"
#include "timer.hpp"

using namespace osh;

static void add_dye(Mesh* mesh) {
  auto dye_w = Write<Real>(mesh->nverts());
  auto coords = mesh->coords();
  auto dye_fun = LAMBDA(LO vert) {
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
  mesh->add_tag(VERT, "dye", 1, OSH_LINEAR_INTERP, Reals(dye_w));
}

static void add_pointwise(Mesh* mesh) {
  auto dim = mesh->dim();
  auto ecoords =
      average_field(mesh, dim, LOs(mesh->nelems(), 0, 1), dim, mesh->coords());
  auto pw_w = Write<Real>(mesh->nelems());
  auto pw_fun = LAMBDA(LO elem) { pw_w[elem] = ecoords[elem * dim]; };
  parallel_for(mesh->nelems(), pw_fun);
  mesh->add_tag(dim, "pointwise", 1, OSH_POINTWISE, Reals(pw_w));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  constexpr Int dim = 3;
  Mesh mesh;
  if (world->rank() == 0) {
    auto nx = 10;
    build_box(&mesh, lib, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
    classify_by_angles(&mesh, PI / 4);
    mesh.reorder();
    mesh.reset_globals();
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_parting(OSH_GHOSTED);
  auto size = find_identity_size(&mesh);
  mesh.add_tag(VERT, "size", 1, OSH_LINEAR_INTERP, size);
  add_dye(&mesh);
  mesh.add_tag(mesh.dim(), "mass", 1, OSH_CONSERVE,
               measure_elements_real(&mesh));
  add_pointwise(&mesh);
  auto mid = zero_vector<dim>();
  mid[0] = mid[1] = .5;
  Now t0 = now();
  for (Int i = 0; i < 8; ++i) {
    auto coords = mesh.coords();
    Write<Real> warp_w(mesh.nverts() * dim);
    auto warp_fun = LAMBDA(LO vert) {
      auto x0 = get_vector<dim>(coords, vert);
      auto x1 = zero_vector<dim>();
      x1[0] = x0[0];
      x1[1] = x0[1];
      auto x2 = x1 - mid;
      auto polar_a = atan2(x2[1], x2[0]);
      auto polar_r = norm(x2);
      Real rot_a = 0;
      if (polar_r < 0.5) {
        rot_a = (PI / 8) * (2.0 * (0.5 - polar_r));
        if (i >= 4) rot_a = -rot_a;
      }
      auto dest_a = polar_a + rot_a;
      auto dst = x0;
      dst[0] = cos(dest_a) * polar_r;
      dst[1] = sin(dest_a) * polar_r;
      dst = dst + mid;
      auto w = dst - x0;
      set_vector<dim>(warp_w, vert, w);
    };
    parallel_for(mesh.nverts(), warp_fun);
    mesh.add_tag(VERT, "warp", dim, OSH_LINEAR_INTERP, Reals(warp_w));
    while (warp_to_limit(&mesh, 0.20)) {
      adapt(&mesh, 0.30, 0.30, 1.0 / 2.0, 3.0 / 2.0, 4, 0);
    }
  }
  Now t1 = now();
  mesh.set_parting(OSH_ELEM_BASED);
  CHECK(are_close(1.0,
                  sum(mesh.comm(), mesh.get_array<Real>(mesh.dim(), "mass"))));
  if (mesh.comm()->rank() == 0) {
    std::cout << "test took " << (t1 - t0) << " seconds\n";
  }
  bool ok = check_regression("gold_warp", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}
