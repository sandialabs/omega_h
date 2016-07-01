#include "internal.hpp"

using namespace osh;

struct Case {
  virtual const char* file_name() const = 0;
  virtual std::vector<I32> objects() const = 0;
  virtual Int time_steps() const = 0;
  virtual Reals motion(Mesh* m, Int step, Int object, LOs object_verts);
};

struct TranslateSphere : public Case {
  virtual const char* file_name() const {
    return "disk_in_square.msh";
  }
  virtual std::vector<I32> objects() const {
    return std::vector<I32>({72});
  }
  virtual Int time_steps() const {
    return 12;
  }
  virtual Reals motion(Mesh* m, Int step, Int object, LOs object_verts) {
    auto out = Write<Real>(object_verts.size() * 3);
    auto f = LAMBDA(LO ov) {
      auto v = object_verts[ov];
      set_vector<3>(out, ov, vector_3(0.02, 0, 0));
    };
    return out;
  }
};

static void run_case(Library const& lib, Case const& c) {
  auto world = lib.world();
  Mesh mesh;
  if (world->rank() == 0) {
    gmsh::read(c.file_name(), lib, &mesh);
  }
  mesh.set_comm(world);
  mesh.balance();
  mesh.set_partition(GHOSTED);
  auto size = find_identity_size(&mesh);
  mesh.add_tag(VERT, "size", 1, OSH_LINEAR_INTERP, size);
  for (Int i = 0; i < c.time_steps(); ++i) {
  }
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
  mesh.set_partition(GHOSTED);
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
      auto x0 = get_vec<dim>(coords, vert);
      auto x1 = zero_vector<dim>();
      x1[0] = x0[0]; x1[1] = x0[1];
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
    while (warp_to_limit(&mesh, 0.30)) {
      adapt(&mesh, 0.30, 0.40, 1.0 / 2.0, 3.0 / 2.0, 4, 0);
    }
  }
  Now t1 = now();
  mesh.set_partition(ELEMENT_BASED);
  CHECK(are_close(1.0, sum(mesh.comm(),
          mesh.get_array<Real>(mesh.dim(), "mass"))));
  if (mesh.comm()->rank() == 0) {
    std::cout << "test took " << (t1 - t0) << " seconds\n";
  }
  bool ok = check_regression("gold_warp", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}

