#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_refine.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_reduce.hpp>

#include <cmath>

using namespace Omega_h;

int main(int argc, char** argv) {
  Library lib(&argc, &argv);
  CommPtr world = lib.world();
  Mesh mesh = build_box(world, OMEGA_H_SIMPLEX, 1., 2., 3., 2, 2, 2);
  AdaptOpts opts(&mesh);
  mesh.add_tag<Real>(VERT, "metric", 1);


  Reals coords  = mesh.coords();
  Write<Real> flux_w(mesh.nfaces());
  {
    auto f = OMEGA_H_LAMBDA(LO v) {
      Vector<3> x = get_vector<3>(coords, v);
      Vector<3> y = {{2*x[0], 5*x[1], x[2]}};
      for (int i=0; i<3; ++i) flux_w[v] = y[0];
    };
    parallel_for(mesh.nfaces(), f);
  }

  const std::string name = "magnetic face flux";
  Reals flux(flux_w);
  mesh.add_tag<Real>(FACE, name, 1, flux);

  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.3)));
  while (coarsen_by_size(&mesh, opts))
    ;

  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.3)));
  while (refine_by_size(&mesh, opts))
    ;

  Read<Real> flux_r=mesh.get_array<Real>(FACE, name); 

  Write<Real> OK(1.,0);
  {
    auto f = OMEGA_H_LAMBDA(LO v) {
      Real x = flux_r[v];
      Vector<3> y = {{2*x, 5*x, x}};
      if (x != y[0]) OK[0] += 1;
      if (x != y[1]) OK[0] += 1;
      if (x != y[2]) OK[0] += 1;
    };
    parallel_for(mesh.nfaces(), f);
  }
  const bool ok = 0.==Reals(OK)[0];

  mesh.ask_qualities();
  if (ok) return 2;
  return 0;
}
