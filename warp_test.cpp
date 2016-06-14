#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  constexpr Int dim = 2;
  Mesh mesh;
  if (world->rank() == 0) {
    auto nx = 8;
    build_box(&mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_partition(GHOSTED);
  auto size = find_identity_size(&mesh);
  mesh.add_tag(VERT, "size", 1, OSH_LINEAR_INTERP, size);
  mesh.add_tag<Real>(VERT, "warp", dim, OSH_LINEAR_INTERP);
  vtk::FullWriter writer(&mesh, "out");
  auto mid = zero_vector<dim>();
  mid[0] = mid[1] = .5;
  for (Int i = 0; i < 1; ++i) {
    auto coords = mesh.coords();
    Write<Real> warp_w(mesh.nverts() * 3);
    auto f = LAMBDA(LO vert) {
      auto x0 = get_vec<dim>(coords, vert);
      auto x1 = zero_vector<dim>();
      x1[0] = x0[0]; x1[1] = x0[1];
      auto polar_a = atan2(x1[1], x1[0]);
      auto polar_r = norm(x1);
      Real rot_a = 0;
      if (polar_r < 0.5) {
        rot_a = (PI / 4) * (2.0 * (0.5 - polar_r));
      }
      auto dest_a = polar_a + rot_a;
      auto dst = zero_vector<dim>();
      dst[0] = cos(dest_a) * polar_r;
      dst[1] = sin(dest_a) * polar_r;
      auto w = dst - x0;
      set_vec<dim>(warp_w, vert, w);
    };
    parallel_for(mesh.nverts(), f);
    mesh.set_tag(VERT, "warp", Reals(warp_w));
    while (warp_to_limit(&mesh, 0.37)) {
      adapt(&mesh, 0.47, 0.37, 2.0 / 3.0, 4.0 / 3.0, 4);
      writer.write();
    }
  }
}

