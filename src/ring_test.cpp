#include "Omega_h.hpp"
#include "Omega_h_math.hpp"
#include "internal.hpp"
#include "laplace.hpp"
#include "loop.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  CHECK(argc == 2);
  auto world = lib.world();
  auto self = lib.self();
  Mesh mesh(&lib);
  gmsh::read(argv[1], &mesh);
  auto ids = std::vector<I32>({6, 7, 8, 9});
  auto verts_are_bdry =
      mark_class_closures(&mesh, VERT, std::vector<Int>(ids.size(), 1), ids);
  auto bv2v = collect_marked(verts_are_bdry);
  auto initial_w = Write<Real>(mesh.nverts(), 0.0);
  auto coords = mesh.coords();
  auto f = LAMBDA(LO bv) {
    auto v = bv2v[bv];
    auto x = get_vector<2>(coords, v);
    auto theta = atan2(x[1], x[0]);
    initial_w[v] = 4.0 * sin(5.0 * theta);
  };
  parallel_for(bv2v.size(), f);
  auto initial = Reals(initial_w);
  auto solution = solve_laplacian(&mesh, initial, 1, 1e-2);
  mesh.add_tag(
      VERT, "solution", 1, OMEGA_H_LINEAR_INTERP, OMEGA_H_DO_OUTPUT, solution);
  bool ok = check_regression("gold_ring", &mesh, 0.0, 0.0);
  if (!ok) return 2;
  return 0;
}
