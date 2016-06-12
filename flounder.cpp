#include "internal.hpp"

using namespace osh;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  Mesh mesh;
  build_box(&mesh, 1, 1, 0, 1, 1, 0);
  classify_by_angles(&mesh, PI / 4);
  mesh.keep_canonical_globals(false);
  mesh.add_tag<Real>(VERT, "size", 1, OSH_LINEAR_INTERP);
  Now t0 = now();
  Int i = 0;
  do {
    Write<Real> size(mesh.nverts());
    auto coords = mesh.coords();
    auto f = LAMBDA(LO v) {
      auto x = get_vec<2>(coords, v);
      auto s = cos(x[0] * 8.0 * M_PI) / 4.0 + 1.0 / 2.0;
      auto d = fabs(x[1] - s);
      auto area = 1e-7 + d * 1e-5;
      auto length = sqrt(2 * area);
      size[v] = length;
    };
    parallel_for(mesh.nverts(), f);
    mesh.set_tag(VERT, "size", Reals(size));
    ++i;
  } while(refine_by_size(&mesh, 0.3));
  Now t1 = now();
  std::cout << "refinement time " << (t1-t0) << " seconds ";
  std::cout << mesh.nelems() << " final triangles\n";
  std::cout << "took " << i << " iterations\n";
}
