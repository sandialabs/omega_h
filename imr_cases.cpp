#include "internal.hpp"

using namespace osh;

struct Case {
  virtual const char* file_name() const = 0;
  virtual std::vector<I32> objects() const = 0;
  virtual Int time_steps() const = 0;
  virtual Reals motion(Mesh* m, Int step, Int object, LOs object_verts) const;
};

struct TranslateBall : public Case {
  virtual const char* file_name() const {
    return "ball_in_cube.msh";
  }
  virtual std::vector<I32> objects() const {
    return std::vector<I32>({72});
  }
  virtual Int time_steps() const {
    return 12;
  }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs object_verts) const {
    (void) m;
    (void) step;
    (void) object;
    auto out = Write<Real>(object_verts.size() * 3);
    auto f = LAMBDA(LO ov) {
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
  mesh.add_tag<Real>(VERT, "warp", mesh.dim(), OSH_LINEAR_INTERP);
  vtk::Writer writer(&mesh, "out", mesh.dim());
  for (Int step = 0; step < c.time_steps(); ++step) {
    auto objs = c.objects();
    auto motion_w = Write<Real>(mesh.nverts() * mesh.dim(), 0.0);
    for (auto obj : objs) {
      auto verts_on_obj = mark_class_closure(
          &mesh, osh::VERT, mesh.dim(), obj);
      auto ov2v = collect_marked(verts_on_obj);
      auto obj_motion = c.motion(&mesh, step, obj, ov2v);
      map_into(obj_motion, ov2v, motion_w, mesh.dim());
    }
    auto motion = Reals(motion_w);
    motion = solve_laplacian(&mesh, motion, mesh.dim(), 1e-3);
    mesh.set_tag(VERT, "warp", motion);
    while (warp_to_limit(&mesh, 0.20)) {
      adapt(&mesh, 0.20, 0.30, 1.0 / 2.0, 3.0 / 2.0, 4, 2);
      writer.write();
    }
  }
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  CHECK(argc == 2);
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
  std::string name = argv[1];
  if (name == "translate_ball") run_case(lib, TranslateBall());
  else osh_fail("unknown case \"%s\"", argv[1]);
}

