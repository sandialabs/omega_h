#include "internal.hpp"

using namespace osh;

struct Case {
  virtual ~Case();
  virtual const char* file_name() const = 0;
  virtual std::vector<I32> objects() const = 0;
  virtual Int time_steps() const = 0;
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const = 0;
};

Case::~Case() {
}

struct TranslateBall : public Case {
  ~TranslateBall();
  virtual const char* file_name() const override {
    return "ball_in_cube.msh";
  }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({72});
  }
  virtual Int time_steps() const override {
    return 12;
  }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void) m;
    (void) step;
    (void) object;
    auto out = Write<Real>(ov2v.size() * 3);
    auto f = LAMBDA(LO ov) {
      set_vector<3>(out, ov, vector_3(0.02, 0, 0));
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

TranslateBall::~TranslateBall() {
}

struct RotateBall : public Case {
  ~RotateBall();
  virtual const char* file_name() const override {
    return "ball_in_cube.msh";
  }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({72});
  }
  virtual Int time_steps() const override {
    return 16;
  }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void) step;
    (void) object;
    auto coords = m->coords();
    auto out = Write<Real>(ov2v.size() * 3);
    auto rot = rotate(PI / 16, vector_3(0,0,1));
    auto f = LAMBDA(LO ov) {
      auto v = ov2v[ov];
      auto x = get_vector<3>(coords, v);
      auto mid = vector_3(.5, .5, 0);
      x = x - mid;
      auto x2 = rot * x;
      auto w = x2 - x;
      set_vector<3>(out, ov, w);
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

RotateBall::~RotateBall() {
}

struct CollideBalls : public Case {
  ~CollideBalls();
  virtual const char* file_name() const override {
    return "balls_in_box.msh";
  }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({72,110});
  }
  virtual Int time_steps() const override {
    return 12;
  }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void) m;
    (void) step;
    auto out = Write<Real>(ov2v.size() * 3);
    auto f = LAMBDA(LO ov) {
      if (object == 72) {
        set_vector<3>(out, ov, vector_3(0, 0, 0.02));
      } else {
        set_vector<3>(out, ov, vector_3(0, 0,-0.02));
      }
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

CollideBalls::~CollideBalls() {
}

struct CylinderTube : public Case {
  ~CylinderTube();
  virtual const char* file_name() const override {
    return "cylinder_thru_tube.msh";
  }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({73});
  }
  virtual Int time_steps() const override {
    return 12;
  }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void) m;
    (void) step;
    (void) object;
    auto out = Write<Real>(ov2v.size() * 3);
    auto f = LAMBDA(LO ov) {
      set_vector<3>(out, ov, vector_3(0, 0, 0.02));
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

CylinderTube::~CylinderTube() {
}

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
  vtk::Writer writer(&mesh, "out", mesh.dim());
  for (Int step = 0; step < c.time_steps(); ++step) {
    mesh.set_partition(GHOSTED);
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
    mesh.add_tag(VERT, "warp", mesh.dim(), OSH_LINEAR_INTERP, motion);
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
  else if (name == "rotate_ball") run_case(lib, RotateBall());
  else if (name == "collide_balls") run_case(lib, CollideBalls());
  else if (name == "cylinder_thru_tube") run_case(lib, CylinderTube());
  else osh_fail("unknown case \"%s\"", argv[1]);
}

