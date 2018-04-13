#include <Omega_h_amr_topology.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_mark.hpp>
#include <Omega_h_mesh.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  mesh = Omega_h::build_box(world, OMEGA_H_HYPERCUBE, 2.0, 1.0, 0.0, 2, 1, 0);
  Omega_h::Write<Omega_h::Byte> elem_mark_w(2, 0);
  elem_mark_w.set(0, 1);
  Omega_h::Read<Omega_h::Byte> elem_mark;
  elem_mark = elem_mark_w;
  Omega_h::mark_amr(&mesh, elem_mark);
  Omega_h::LOs empty;
  Omega_h::LOs p2mv0({0, 2, 9, 7});
  Omega_h::LOs p2mv1({1, 6, 8, 4});
  Omega_h::LOs p2mv2({5});
  Omega_h::Few<Omega_h::LOs, 4> p2mv({p2mv0, p2mv1, p2mv2, empty});
  Omega_h::LOs p2mds0({0, 4, 2, 1});
  Omega_h::LOs p2mds1({0, 2, 3, 1});
  Omega_h::LOs p2mds2({0});
  Omega_h::Few<Omega_h::LOs, 4> p2mds({p2mds0, p2mds1, p2mds2, empty});
  Omega_h::Few<Omega_h::LOs, 4> mds2p;
  for (Omega_h::Int dim = 0; dim <= mesh.dim(); ++dim) {
    if (p2mds[dim].exists()) {
      mds2p[dim] = Omega_h::invert_injective_map(p2mds[dim], mesh.nents(dim));
    }
  }
  auto counts = Omega_h::count_amr(&mesh);
  auto num_child_edges = counts[1];
  auto edge_verts =
      Omega_h::get_amr_topology(&mesh, 1, num_child_edges, p2mds, mds2p, p2mv);
  Omega_h::LOs truth(
      {0, 1, 1, 2, 9, 6, 6, 2, 9, 8, 8, 7, 7, 4, 4, 0, 1, 5, 6, 5, 8, 5, 4, 5});
  OMEGA_H_CHECK(truth.size() == 24);
  OMEGA_H_CHECK(edge_verts == truth);
  return 0;
}
