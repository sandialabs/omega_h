#include "Omega_h_adapt.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_laplace.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_vector.hpp"

#include <iostream>
#include <set>

using namespace Omega_h;

struct Case {
  virtual ~Case();
  virtual const char* file_name() const = 0;
  virtual std::vector<I32> objects() const = 0;
  virtual Int time_steps() const = 0;
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const = 0;
};

Case::~Case() {}

struct TranslateBall : public Case {
  ~TranslateBall() override;
  virtual const char* file_name() const override { return "ball_in_cube.msh"; }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({72});
  }
  virtual Int time_steps() const override { return 12; }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void)m;
    (void)step;
    (void)object;
    return static_motion(ov2v);
  }
  static Reals static_motion(LOs ov2v) {
    auto out = Write<Real>(ov2v.size() * 3);
    auto f = OMEGA_H_LAMBDA(LO ov) {
      set_vector<3>(out, ov, vector_3(0.02, 0, 0));
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

TranslateBall::~TranslateBall() {}

struct RotateBall : public Case {
  ~RotateBall() override;
  virtual const char* file_name() const override { return "ball_in_cube.msh"; }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({72});
  }
  virtual Int time_steps() const override { return 16; }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void)step;
    (void)object;
    return static_motion(m, ov2v);
  }
  static Reals static_motion(Mesh* m, LOs ov2v) {
    auto coords = m->coords();
    auto out = Write<Real>(ov2v.size() * 3);
    auto rot = rotate(PI / 16, vector_3(0, 0, 1));
    auto f = OMEGA_H_LAMBDA(LO ov) {
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

RotateBall::~RotateBall() {}

struct CollideBalls : public Case {
  ~CollideBalls() override;
  virtual const char* file_name() const override { return "balls_in_box.msh"; }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({72, 110});
  }
  virtual Int time_steps() const override { return 12; }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void)m;
    (void)step;
    return static_motion(object, ov2v);
  }
  static Reals static_motion(I32 object, LOs ov2v) {
    auto out = Write<Real>(ov2v.size() * 3);
    auto f = OMEGA_H_LAMBDA(LO ov) {
      if (object == 72) {
        set_vector<3>(out, ov, vector_3(0, 0, 0.02));
      } else {
        set_vector<3>(out, ov, vector_3(0, 0, -0.02));
      }
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

CollideBalls::~CollideBalls() {}

struct CylinderTube : public Case {
  ~CylinderTube() override;
  virtual const char* file_name() const override {
    return "cylinder_thru_tube.msh";
  }
  virtual std::vector<I32> objects() const override {
    return std::vector<I32>({73});
  }
  virtual Int time_steps() const override { return 12; }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void)m;
    (void)step;
    (void)object;
    return static_motion(ov2v);
  }
  static Reals static_motion(LOs ov2v) {
    auto out = Write<Real>(ov2v.size() * 3);
    auto f = OMEGA_H_LAMBDA(LO ov) {
      set_vector<3>(out, ov, vector_3(0, 0, 0.02));
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

CylinderTube::~CylinderTube() {}

struct TwinRotor : public Case {
  std::set<I32> assembly0;
  std::set<I32> assembly1;
  TwinRotor() : assembly0({66, 98, 126}), assembly1({254, 253, 252}) {}
  ~TwinRotor() override;
  virtual const char* file_name() const override { return "twin_rotor.msh"; }
  virtual std::vector<I32> objects() const override {
    std::vector<I32> out;
    out.insert(out.end(), assembly0.begin(), assembly0.end());
    out.insert(out.end(), assembly1.begin(), assembly1.end());
    return out;
  }
  virtual Int time_steps() const override { return 16; }
  virtual Reals motion(Mesh* m, Int step, I32 object, LOs ov2v) const override {
    (void)step;
    Vector<3> center;
    Real dir;
    if (assembly0.count(object)) {
      center = vector_3(-.25, 0, 0);
      dir = 1.0;
    } else if (assembly1.count(object)) {
      center = vector_3(.25, 0, 0);
      dir = -1.0;
    } else {
      Omega_h_fail("object %d not in either assembly\n", object);
    }
    return static_motion(m, ov2v, center, dir);
  }
  static Reals static_motion(Mesh* m, LOs ov2v, Vector<3> center, Real dir) {
    auto coords = m->coords();
    auto out = Write<Real>(ov2v.size() * 3);
    auto rm = rotate(dir * PI / 32, vector_3(0, 0, 1));
    auto f = OMEGA_H_LAMBDA(LO ov) {
      auto v = ov2v[ov];
      auto x = get_vector<3>(coords, v);
      set_vector(out, ov, ((rm * (x - center)) + center) - x);
    };
    parallel_for(ov2v.size(), f);
    return out;
  }
};

TwinRotor::~TwinRotor() {}

static void run_case(Library* lib, Case const& c, Int niters) {
  if (niters == -1) {
    niters = c.time_steps();
  } else {
    OMEGA_H_CHECK(niters >= 0);
    if (niters > c.time_steps()) {
      std::cerr << "warning: requesting " << niters
                << " time steps but the case is designed for " << c.time_steps()
                << '\n';
    }
  }
  auto world = lib->world();
  auto mesh = gmsh::read(c.file_name(), world);
  mesh.set_parting(OMEGA_H_GHOSTED);
  {
    auto metrics = get_implied_isos(&mesh);
    mesh.add_tag(VERT, "metric", 1, metrics);
  }
  vtk::Writer writer("out", &mesh);
  writer.write();
  Now t0 = now();
  for (Int step = 0; step < niters; ++step) {
    mesh.set_parting(OMEGA_H_GHOSTED);
    auto objs = c.objects();
    auto motion_w = Write<Real>(mesh.nverts() * mesh.dim(), 0.0);
    for (auto obj : objs) {
      auto verts_on_obj =
          mark_class_closure(&mesh, Omega_h::VERT, mesh.dim(), obj);
      auto ov2v = collect_marked(verts_on_obj);
      auto obj_motion = c.motion(&mesh, step, obj, ov2v);
      map_into(obj_motion, ov2v, motion_w, mesh.dim());
    }
    auto motion = Reals(motion_w);
    motion = solve_laplacian(&mesh, motion, mesh.dim(), 1e-2);
    mesh.add_tag(VERT, "warp", mesh.dim(), motion);
    // auto metrics = mesh.get_array<Real>(VERT, "metric");
    // auto lengths = lengths_from_isos(metrics);
    // lengths = solve_laplacian(&mesh, lengths, 1, 1e-2);
    // metrics = isos_from_lengths(lengths);
    // mesh.set_tag(VERT, "metric", metrics);
    auto opts = AdaptOpts(&mesh);
    opts.max_length_allowed = opts.max_length_desired * 2.0;
    int warp_step = 0;
    while (warp_to_limit(&mesh, opts)) {
      if (world->rank() == 0) {
        std::cout << "WARP STEP " << warp_step << " OF TIME STEP " << step
                  << '\n';
      }
      adapt(&mesh, opts);
      ++warp_step;
    }
    writer.write();
  }
  Now t1 = now();
  if (world->rank() == 0) {
    std::cout << "case took " << (t1 - t0) << " seconds to run\n";
  }
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  std::string name;
  Int niters = -1;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--niters") {
      if (i == argc - 1) Omega_h_fail("--niters needs an argument\n");
      ++i;
      niters = atoi(argv[i]);
    } else {
      name = arg;
    }
  }
  if (name == "translate_ball")
    run_case(&lib, TranslateBall(), niters);
  else if (name == "rotate_ball")
    run_case(&lib, RotateBall(), niters);
  else if (name == "collide_balls")
    run_case(&lib, CollideBalls(), niters);
  else if (name == "cylinder_thru_tube")
    run_case(&lib, CylinderTube(), niters);
  else if (name == "twin_rotor")
    run_case(&lib, TwinRotor(), niters);
  else
    Omega_h_fail("unknown case \"%s\"", argv[1]);
}
