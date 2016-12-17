#include "Omega_h.hpp"
#include "Omega_h_math.hpp"
#include "simplices.hpp"
#include "loop.hpp"
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
  auto f = LAMBDA(LO v) {
    auto p = get_vector<dim>(coords, v);
    auto x = p[0];
    auto y = p[1];
    auto xy = vector_2(x, y);
    auto radius = norm(xy);
    auto t = atan2(y, x);
    auto h = vector_3(h0 + 2 * (0.1 - h0) * fabs(radius - 0.5), 0.1, 0.1);
    auto rotation = rotate(t, vector_3(0,0,1));
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
  Omega_h_fail("no UGAWG metric named %s\n", name.c_str());
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto opts = AdaptOpts(dim);
  std::string filename;
  std::string metric_name;
  bool should_help = false;
  for (Int i = 1; i < argc; ++i) {
    if (std::string("-m") == argv[i]) {
      if (i == argc - 1) {
        std::cout << "-m needs an argument\n";
        should_help = true;
        break;
      }
      metric_name = argv[++i];
    } else if (std::string("-h") == argv[i]) {
      if (i >= argc - 2) {
        std::cout << "-h needs two arguments\n";
        should_help = true;
        break;
      }
      opts.length_histogram_min = atof(argv[++i]);
      opts.length_histogram_max = atof(argv[++i]);
    } else if (std::string("-l") == argv[i]) {
      if (i >= argc - 2) {
        std::cout << "-l needs two arguments\n";
        should_help = true;
        break;
      }
      opts.min_length_desired = atof(argv[++i]);
      opts.max_length_desired = atof(argv[++i]);
    } else if (std::string("-q") == argv[i]) {
      if (i >= argc - 1) {
        std::cout << "-q needs two arguments\n";
        should_help = true;
        break;
      }
      opts.min_quality_desired = atof(argv[++i]);
    } else {
      filename = argv[i];
    }
  }
  if (metric_name.empty() || filename.empty()) should_help = true;
  if (should_help) {
    std::cout << "usage: " << argv[0] << " [options] input.mesh[b]\n";
    std::cout << "options:\n";
    std::cout << "  -m <metric-name>       (REQUIRED) one of:\n";
    std::cout << "                            cube-linear\n";
    std::cout << "                            cube-cylinder-shock\n";
    std::cout << "                            cube-cylinder-layer\n";
    std::cout << "  -h <length-histogram-min> <length-histogram-min>\n";
    std::cout << "  -l <min-desired-length> <max-desired-length>\n";
    std::cout << "  -q <min-desired-quality>\n";
    return -1;
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::meshb::read(&mesh, filename.c_str());
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
