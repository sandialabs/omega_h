#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  constexpr Int dim = 2;
  Mesh mesh;
  if (world->rank() == 0) {
    auto nx = 12;
    build_box(&mesh, 1, 1, 1, nx, nx, (dim == 3) ? nx : 0);
    classify_by_angles(&mesh, PI / 4);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_partition(GHOSTED);
  auto size = find_identity_size(&mesh);
  mesh.add_tag(VERT, "size", 1, OSH_LINEAR_INTERP, size);
  vtk::FullWriter writer(&mesh, "out");
  auto mid = zero_vector<dim>();
  mid[0] = mid[1] = .5;
  for (Int i = 0; i < 4; ++i) {
    if (world->rank() == 0) std::cout << "rotation step " << i << '\n';
    auto coords = mesh.coords();
    Write<Real> warp_w(mesh.nverts() * dim);
    auto f = LAMBDA(LO vert) {
      auto x0 = get_vec<dim>(coords, vert);
      auto x1 = zero_vector<dim>();
      x1[0] = x0[0]; x1[1] = x0[1];
      auto x2 = x1 - mid;
      auto polar_a = atan2(x2[1], x2[0]);
      auto polar_r = norm(x2);
      Real rot_a = 0;
      if (polar_r < 0.5) {
        rot_a = (PI / 8) * (2.0 * (0.5 - polar_r));
      }
      auto dest_a = polar_a + rot_a;
      auto dst = x0;
      dst[0] = cos(dest_a) * polar_r;
      dst[1] = sin(dest_a) * polar_r;
      dst = dst + mid;
      auto w = dst - x0;
      set_vec<dim>(warp_w, vert, w);
    };
    parallel_for(mesh.nverts(), f);
    mesh.add_tag(VERT, "warp", dim, OSH_LINEAR_INTERP, Reals(warp_w));
    Int wi = 0;
    while (warp_to_limit(&mesh, 0.37)) {
      if (world->rank() == 0) std::cout << "warp step " << wi++ << '\n';
      adapt(&mesh, 0.37, 0.47, 2.0 / 3.0, 4.0 / 3.0, 4);
      writer.write();
    }
  }
}

