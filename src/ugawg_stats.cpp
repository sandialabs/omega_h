#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_element.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_histogram.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_metric.hpp>

#include <iostream>

using namespace Omega_h;

constexpr Int dim = 3;

static Reals get_cube_linear_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto z = coords[v * dim + (dim - 1)];
    auto h = Vector<dim>();
    for (Int i = 0; i < dim - 1; ++i) h[i] = 0.1;
    h[dim - 1] = 0.001 + 0.198 * std::abs(z - 0.5);
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_cube_cylinder_shock_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  constexpr Real h0 = 0.001;
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto z = p[2];
    auto h = vector_3(0.1, 0.1, h0 + 2 * (0.1 - h0) * std::abs(z - 0.5));
    auto m = diagonal(metric_eigenvalues_from_lengths(h));
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_cube_cylinder_layer_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  constexpr Real h0 = 0.001;
  constexpr Real h_z = 1.0 / 10.0;
  constexpr Real h_t = 1.0 / 10.0;
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = std::atan2(y, x);
    auto h = vector_3(h0 + 2 * (0.1 - h0) * std::abs(radius - 0.5), h_t, h_z);
    auto rotation = rotate(t, vector_3(0, 0, 1));
    auto m = compose_metric(rotation, h);
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_cube_cylinder_quality_layer_metric(Mesh* mesh) {
  auto coords = mesh->coords();
  auto out = Write<Real>(mesh->nverts() * symm_ncomps(dim));
  constexpr Real h0 = 0.001;
  constexpr Real h_z = 1.0 / 10.0;
  auto f = OMEGA_H_LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = std::atan2(y, x);
    auto d = (0.6 - radius) * 10.0;
    Real h_t = (d < 0.0) ? (1.0 / 10.0)
                         : (d * (1.0 / 40.0) + (1.0 - d) * (1.0 / 10.0));
    auto h = vector_3(h0 + 2 * (0.1 - h0) * std::abs(radius - 0.5), h_t, h_z);
    auto rotation = rotate(t, vector_3(0, 0, 1));
    auto m = compose_metric(rotation, h);
    set_symm(out, v, m);
  };
  parallel_for(mesh->nverts(), f);
  return out;
}

static Reals get_metric(Mesh* mesh, std::string const& name) {
  if (name == "implied") {
    return get_implied_metrics(mesh);
  }
  if (name == "cube-linear") {
    return get_cube_linear_metric(mesh);
  }
  if (name == "cube-cylinder-linear") {
    return get_cube_cylinder_shock_metric(mesh);
  }
  if (name == "cube-cylinder-polar-1") {
    return get_cube_cylinder_layer_metric(mesh);
  }
  if (name == "cube-cylinder-polar-2") {
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
      "cube-linear,cube-cylinder-linear,cube-cylinder-polar-1,cube-cylinder-"
      "polar-2");
  mflag.add_arg<std::string>("metric");
  auto& hflag = cmdline.add_flag("-h", "domain of length histogram");
  hflag.add_arg<double>("length-histogram-min");
  hflag.add_arg<double>("length-histogram-max");
  auto& nflag = cmdline.add_flag("-n", "numbers of histogram bins");
  nflag.add_arg<int>("quality-histogram-bins");
  nflag.add_arg<int>("length-histogram-bins");
  auto& fflag = cmdline.add_flag("-f", "filenames for histograms");
  fflag.add_arg<std::string>("quality-histogram-file");
  fflag.add_arg<std::string>("length-histogram-file");
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
  if (cmdline.parsed("-n")) {
    opts.nquality_histogram_bins =
        cmdline.get<int>("-n", "quality-histogram-bins");
    opts.nlength_histogram_bins =
        cmdline.get<int>("-n", "length-histogram-bins");
  }
  if (cmdline.parsed("-l")) {
    opts.min_length_desired = cmdline.get<double>("-l", "min-desired-length");
    opts.max_length_desired = cmdline.get<double>("-l", "max-desired-length");
  }
  if (cmdline.parsed("-q")) {
    opts.min_quality_desired = cmdline.get<double>("-q", "min-desired-quality");
  }
  for (Int i = 0; i <= dim; ++i) {
    std::cout << "mesh has " << mesh.nents(i) << ' '
              << topological_plural_name(mesh.family(), i) << '\n';
  }
  auto metrics = get_metric(&mesh, metric_name);
  mesh.add_tag(VERT, "metric", symm_ncomps(dim), metrics);
  print_adapt_status(&mesh, opts);
  auto qh = get_histogram(&mesh, mesh.dim(), opts.nquality_histogram_bins, 0.0,
      1.0, mesh.ask_qualities());
  auto lh = get_histogram(&mesh, EDGE, opts.nlength_histogram_bins,
      opts.length_histogram_min, opts.length_histogram_max, mesh.ask_lengths());
  if (cmdline.parsed("-f")) {
    auto qf = cmdline.get<std::string>("-f", "quality-histogram-file");
    auto lf = cmdline.get<std::string>("-f", "length-histogram-file");
    render_histogram_matplotlib(qh, qf);
    render_histogram_matplotlib(lh, lf);
  } else {
    print_histogram(qh, "quality");
    print_histogram(lh, "length");
  }
  return 0;
}
