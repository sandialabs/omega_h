#include "Omega_h.hpp"
#include "Omega_h_cmdline.hpp"
#include "Omega_h_math.hpp"
#include "loop.hpp"
#include "simplices.hpp"
#include "space.hpp"

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 3;

static Reals get_cube_linear_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_dofs(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.1;
    h[dim - 1] = 0.001 + 0.198 * fabs(z - 0.5);
    auto m = diagonal(metric_eigenvalues(h));
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_cube_cylinder_shock_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_dofs(dim));
  constexpr Real h0 = 0.001;
  auto f = LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto z = p[2];
    auto h = vector_3(0.1, 0.1, h0 + 2 * (0.1 - h0) * fabs(z - 0.5));
    auto m = diagonal(metric_eigenvalues(h));
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_cube_cylinder_layer_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_dofs(dim));
  constexpr Real h0 = 0.001;
  constexpr Real h_z = 1.0 / 10.0;
  constexpr Real h_t = 1.0 / 10.0;
  auto f = LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = atan2(y, x);
    auto h = vector_3(h0 + 2 * (0.1 - h0) * fabs(radius - 0.5), h_t, h_z);
    auto rotation = rotate(t, vector_3(0, 0, 1));
    auto m = compose_metric(rotation, h);
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_cube_cylinder_quality_layer_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_dofs(dim));
  constexpr Real h0 = 0.001;
  constexpr Real h_z = 1.0 / 10.0;
  auto f = LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = atan2(y, x);
    auto d = (0.6 - radius) * 10.0;
    Real h_t = (d < 0.0) ? (1.0 / 10.0)
                         : (d * (1.0 / 40.0) + (1.0 - d) * (1.0 / 10.0));
    auto h = vector_3(h0 + 2 * (0.1 - h0) * fabs(radius - 0.5), h_t, h_z);
    auto rotation = rotate(t, vector_3(0, 0, 1));
    auto m = compose_metric(rotation, h);
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_metric(Mesh* mesh, std::string const& name) {
  if (name == "cube-linear") {
    return get_cube_linear_metric(mesh);
  }
  if (name == "cube-cylinder-shock") {
    return get_cube_cylinder_shock_metric(mesh);
  }
  if (name == "cube-cylinder-layer") {
    return get_cube_cylinder_layer_metric(mesh);
  }
  if (name == "cube-cylinder-quality-layer") {
    return get_cube_cylinder_quality_layer_metric(mesh);
  }
  Omega_h_fail("no UGAWG metric named %s\n", name.c_str());
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("input.mesh[b]");
  auto& mflag = cmdline.add_flag("-m",
      "REQUIRED!\n    one of "
      "cube-linear,cube-cylinder-shock,cube-cylinder-layer");
  mflag.add_arg<std::string>("metric");
  auto& hflag = cmdline.add_flag("-h", "domain of length histogram");
  hflag.add_arg<double>("length-histogram-min");
  hflag.add_arg<double>("length-histogram-max");
  auto& lflag = cmdline.add_flag("-l", "range of desired lengths");
  lflag.add_arg<double>("min-desired-length");
  lflag.add_arg<double>("max-desired-length");
  auto& qflag = cmdline.add_flag("-q", "quality target");
  qflag.add_arg<double>("min-desired-quality");
  auto world = lib.world();
  if (!cmdline.parse(world, &argc, argv)) {
    std::cout << "parse failed\n";
    cmdline.show_help(world, argv);
    return -1;
  }
  if (!cmdline.check_empty(world, argc, argv)) {
    std::cout << "leftover args\n";
    cmdline.show_help(world, argv);
    return -1;
  }
  if (!cmdline.parsed("-m")) {
    std::cout << "-m wasn't parsed\n";
    cmdline.show_help(world, argv);
    return -1;
  }
  std::string filename = cmdline.get<std::string>("input.mesh[b]");
  std::string metric_name = cmdline.get<std::string>("-m", "metric");
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, filename.c_str());
  auto opts = AdaptOpts(&mesh);
  if (cmdline.parsed("-h")) {
    opts.length_histogram_min =
        cmdline.get<double>("-h", "length-histogram-min");
    opts.length_histogram_max =
        cmdline.get<double>("-h", "length-histogram-max");
  }
  if (cmdline.parsed("-l")) {
    opts.min_length_desired = cmdline.get<double>("-l", "min-desired-length");
    opts.max_length_desired = cmdline.get<double>("-l", "max-desired-length");
  }
  if (cmdline.parsed("-q")) {
    opts.min_quality_desired = cmdline.get<double>("-q", "min-desired-quality");
  }
  for (Int i = 0; i <= dim; ++i) {
    std::cout << "mesh has " << mesh.nents(i) << ' ' << plural_names[i] << '\n';
  }
  auto metrics = get_metric(&mesh, metric_name);
  mesh.add_tag(VERT, "metric", symm_dofs(dim), OMEGA_H_METRIC,
      OMEGA_H_DO_OUTPUT, metrics);
  print_adapt_status(&mesh, opts);
  print_adapt_histograms(&mesh, opts);
  return 0;
}
