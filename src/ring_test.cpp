#include "Omega_h.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_internal.hpp"
#include "Omega_h_laplace.hpp"
#include "Omega_h_map.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  CHECK(argc == 2);
  auto world = lib.world();
  auto self = lib.self();
  Mesh mesh(&lib);
  gmsh::read(argv[1], &mesh);
  auto ids = std::vector<LO>({6, 7, 8, 9});
  auto verts_are_bdry =
      mark_class_closures(&mesh, VERT, 1, ids);
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
  mesh.add_tag(VERT, "solution", 1, solution);
  bool ok = check_regression("gold_ring", &mesh);
  if (!ok) return 2;
  return 0;
}
