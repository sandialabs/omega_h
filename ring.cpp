#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  auto self = lib.self();
  Mesh mesh;
  gmsh::read("ring.msh", &mesh);
  auto ids = std::vector<I32>({6,7,8,9});
  auto verts_are_bdry = mark_class_closures(&mesh, VERT,
      std::vector<Int>(ids.size(), 1), ids);
  auto bv2v = collect_marked(verts_are_bdry);
  auto initial_w = Write<Real>(mesh.nverts(), 0.0);
  auto coords = mesh.coords();
  auto f = LAMBDA(LO bv) {
    auto v = bv2v[bv];
    auto x = get_vec<2>(coords, v);
    auto theta = atan2(x[1], x[0]);
    initial_w[v] = 4.0 * sin(5.0 * theta);
  };
  parallel_for(bv2v.size(), f);
  auto initial = Reals(initial_w);
  auto solution = solve_laplacian(&mesh, initial, 1, 1e-2);
  mesh.add_tag(VERT, "solution", 1, OSH_LINEAR_INTERP, solution);
  vtk::write_parallel("out", &mesh, mesh.dim());
}
