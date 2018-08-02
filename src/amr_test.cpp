#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>

using Omega_h::Bytes;
using Omega_h::Int;
using Omega_h::LOs;

static void check_leaves_before(Omega_h::Mesh* m, Int* lengths) {
  for (Int i = 0; i < m->dim(); ++i) {
    auto leaves = m->ask_leaves(i);
    OMEGA_H_CHECK(leaves == Bytes(lengths[i], 1));
  }
}

static void check_levels_before(Omega_h::Mesh* m, Int* lengths) {
  for (Int i = 0; i < m->dim(); ++i) {
    auto levels = m->ask_levels(i);
    OMEGA_H_CHECK(levels == Bytes(lengths[i], 0));
  }
}

static void check_parents_before(Omega_h::Mesh* m, Int* lengths) {
  for (Int i = 0; i < m->dim(); ++i) {
    auto parents = m->ask_parents(i);
    OMEGA_H_CHECK(parents.parent_idx == LOs(lengths[i], -1));
    OMEGA_H_CHECK(parents.codes == Bytes(lengths[i], 0));
  }
}

static void check_before(Omega_h::Mesh* m, Int* lengths) {
  check_leaves_before(m, lengths);
  check_levels_before(m, lengths);
  check_parents_before(m, lengths);
}

static void check_2D_leaves_after(Omega_h::Mesh* m) {
  auto vtx_leaves = m->ask_leaves(0);
  auto edge_leaves = m->ask_leaves(1);
  auto face_leaves = m->ask_leaves(2);
  auto e = Bytes({0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1});
  OMEGA_H_CHECK(vtx_leaves == Bytes(11, 1));
  OMEGA_H_CHECK(edge_leaves == e);
  OMEGA_H_CHECK(face_leaves == Bytes({0, 1, 1, 1, 1, 1}));
}

static void check_2D_levels_after(Omega_h::Mesh* m) {
  auto vtx_levels = m->ask_levels(0);
  auto edge_levels = m->ask_levels(1);
  auto face_levels = m->ask_levels(2);
  auto v = Bytes({0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0});
  auto e = Bytes({0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0});
  auto f = Bytes({0, 1, 1, 1, 1, 0});
  OMEGA_H_CHECK(vtx_levels == v);
  OMEGA_H_CHECK(edge_levels == e);
  OMEGA_H_CHECK(face_levels == f);
}

static void check_2D_parents_after(Omega_h::Mesh* m) {
  auto vtx_parents = m->ask_parents(0);
  auto edge_parents = m->ask_parents(1);
  auto face_parents = m->ask_parents(2);
  auto vp = LOs({-1, 0, 0, -1, 7, -1, 10, 13, -1, -1, -1});
  auto vc = Bytes({0, 1, 2, 0, 1, 0, 1, 1, 0, 0, 0});
  auto ep =
      LOs({-1, 0, 0, 0, 0, 0, 0, -1, 7, 7, -1, 10, 10, -1, 13, 13, -1, -1, -1});
  auto ec = Bytes({0, 1, 5, 2, 6, 10, 14, 0, 1, 5, 0, 1, 5, 0, 1, 5, 0, 0, 0});
  auto fp = LOs({-1, 0, 0, 0, 0, -1});
  auto fc = Bytes({0, 2, 6, 10, 14, 0});
  OMEGA_H_CHECK(vtx_parents.parent_idx == vp);
  OMEGA_H_CHECK(vtx_parents.codes == vc);
  OMEGA_H_CHECK(edge_parents.parent_idx == ep);
  OMEGA_H_CHECK(edge_parents.codes == ec);
  OMEGA_H_CHECK(face_parents.parent_idx == fp);
  OMEGA_H_CHECK(face_parents.codes == fc);
}

static void check_2D_children_after(Omega_h::Mesh* m) {
  auto edge_vtx_children = m->ask_children(1, 0);
  auto edge_edge_children = m->ask_children(1, 1);
  auto face_vtx_children = m->ask_children(2, 0);
  auto face_edge_children = m->ask_children(2, 1);
  auto face_face_children = m->ask_children(2, 2);
  auto ev_a2ab =
      LOs({0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4});
  auto ev_ab2b = LOs({1, 4, 6, 7});
  auto ev_codes = Bytes({1, 1, 1, 1});
  auto ee_a2ab =
      LOs({0, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 6, 6, 6, 8, 8, 8, 8, 8, 8});
  auto ee_ab2b = LOs({1, 2, 8, 9, 11, 12, 14, 15});
  auto ee_codes = Bytes({1, 5, 1, 5, 1, 5, 1, 5});
  OMEGA_H_CHECK(edge_vtx_children.a2ab == ev_a2ab);
  OMEGA_H_CHECK(edge_vtx_children.ab2b == ev_ab2b);
  OMEGA_H_CHECK(edge_vtx_children.codes == ev_codes);
  OMEGA_H_CHECK(edge_edge_children.a2ab == ee_a2ab);
  OMEGA_H_CHECK(edge_edge_children.ab2b == ee_ab2b);
  OMEGA_H_CHECK(edge_edge_children.codes == ee_codes);
  OMEGA_H_CHECK(face_vtx_children.a2ab == LOs({0, 1, 1, 1, 1, 1, 1}));
  OMEGA_H_CHECK(face_vtx_children.ab2b == LOs({2}));
  OMEGA_H_CHECK(face_vtx_children.codes == Bytes({2}));
  OMEGA_H_CHECK(face_edge_children.a2ab == LOs({0, 4, 4, 4, 4, 4, 4}));
  OMEGA_H_CHECK(face_edge_children.ab2b == LOs({3, 4, 5, 6}));
  OMEGA_H_CHECK(face_edge_children.codes == Bytes({2, 6, 10, 14}));
  OMEGA_H_CHECK(face_face_children.a2ab == LOs({0, 4, 4, 4, 4, 4, 4}));
  OMEGA_H_CHECK(face_face_children.ab2b == LOs({1, 2, 3, 4}));
  OMEGA_H_CHECK(face_face_children.codes == Bytes({2, 6, 10, 14}));
}

static void check_2D_after(Omega_h::Mesh* m) {
  check_2D_leaves_after(m);
  check_2D_levels_after(m);
  check_2D_parents_after(m);
  check_2D_children_after(m);
}

static void test_2D_arrays(Omega_h::Library* lib) {
  auto w = lib->world();
  auto f = OMEGA_H_HYPERCUBE;
  auto m = Omega_h::build_box(w, f, 2.0, 1.0, 0.0, 2, 1, 0);
  Int lengths[3] = {6, 7, 2};
  check_before(&m, lengths);
  auto xfer_opts = Omega_h::TransferOpts();
  Omega_h::amr_refine(&m, Omega_h::Bytes({1, 0}), xfer_opts);
  check_2D_after(&m);
}

static void check_3D_leaves_after(Omega_h::Mesh* m) {
  auto vtx_leaves = m->ask_leaves(0);
  auto edge_leaves = m->ask_leaves(1);
  auto face_leaves = m->ask_leaves(2);
  auto region_leaves = m->ask_leaves(3);
  auto e = Bytes({0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0,
      1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1,
      1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1,
      1, 1, 1});
  auto f = Bytes({0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0,
      1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1,});
  auto r = Bytes({0, 1, 1, 1, 1, 1, 1, 1, 1, 1});
  OMEGA_H_CHECK(vtx_leaves == Bytes(31, 1));
  OMEGA_H_CHECK(edge_leaves == e);
  OMEGA_H_CHECK(face_leaves == f);
  OMEGA_H_CHECK(region_leaves == r);
}

static void check_3D_levels_after(Omega_h::Mesh* m) {
  auto vtx_levels = m->ask_levels(0);
  auto edge_levels = m->ask_levels(1);
  auto face_levels = m->ask_levels(2);
  auto region_levels = m->ask_levels(3);
  auto v = Bytes({0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0,
      0, 1, 1, 0, 0, 0, 0, 1, 1, 0});
  auto e = Bytes({0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0,
      1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1,
      1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1,
      1, 0, 0});
  auto f = Bytes({0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0,
      1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
      0});
  auto r = Bytes({0, 1, 1, 1, 1, 1, 1, 1, 1, 0});
  OMEGA_H_CHECK(vtx_levels == v);
  OMEGA_H_CHECK(edge_levels == e);
  OMEGA_H_CHECK(face_levels == f);
  OMEGA_H_CHECK(region_levels == r);
}

static void check_3D_parents_after(Omega_h::Mesh* m) {
  auto vtx_parents = m->ask_parents(0);
  auto edge_parents = m->ask_parents(1);
  auto face_parents = m->ask_parents(2);
  auto region_parents = m->ask_parents(3);
  auto vp = LOs({-1, 0, 17, 0, 0, -1, 20, 5, -1, 27, -1, 30, 37, 10, -1, 40, 47,
      54, 15, 20, -1, -1, 59, 62, -1, -1, -1, -1, 69, 27, -1});
  auto vc = Bytes({0, 1, 1, 2, 3, 0, 1, 2, 0, 1, 0, 1, 1, 2, 0, 1, 1, 1, 2, 2,
      0, 0, 1, 1, 0, 0, 0, 0, 1, 2, 0});
  auto ep = LOs({-1, 0, 0, 0, 0, 0, 0, 27, 27, 27, 27, 0, 0, 0, 0, 0, 0, -1, 17,
      17, -1, 20, 20, 5, 5, 5, 5, -1, 27, 27, -1, 30, 30, 10, 10, 10, 10, -1,
      37, 37, -1, 40, 40, 15, 15, 15, 15, -1, 47, 47, 20, 20, 20, 20, -1, 54,
      54, -1, -1, -1, 59, 59, -1, 62, 62, -1, -1, -1, -1, -1, 69, 69, -1, -1});
  auto ec = Bytes({0, 1, 5, 2, 6, 10, 14, 2, 6, 10, 14, 3, 7, 11, 15, 19, 23, 0,
      1, 5, 0, 1, 5, 2, 6, 10, 14, 0, 1, 5, 0, 1, 5, 2, 6, 10, 14, 0, 1, 5, 0,
      1, 5, 2, 6, 10, 14, 0, 1, 5, 2, 6, 10, 14, 0, 1, 5, 0, 0, 0, 1, 5, 0, 1,
      5, 0, 0, 0, 0, 0, 1, 5, 0, 0});
  auto fp = LOs({-1, 0, 0, 0, 0, -1, 5, 5, 5, 5, -1, 10, 10, 10, 10, -1, 15, 15,
      15, 15, -1, 20, 20, 20, 20, -1, -1, -1, 27, 27, 27, 27, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, -1, -1, -1});
  auto fc = Bytes({0, 2, 6, 10, 14, 0, 2, 6, 10, 14, 0, 2, 6, 10, 14, 0, 2, 6,
      10, 14, 0, 2, 6, 10, 14, 0, 0, 0, 2, 6, 10, 14, 3, 7, 11, 15, 19, 23, 27,
      31, 35, 39, 43, 47, 0, 0, 0});
  auto rp = LOs({-1, 0, 0, 0, 0, 0, 0, 0, 0, -1});
  auto rc = Bytes({0, 3, 7, 11, 15, 19, 23, 27, 31, 0});
  OMEGA_H_CHECK(vtx_parents.parent_idx == vp);
  OMEGA_H_CHECK(vtx_parents.codes == vc);
  OMEGA_H_CHECK(edge_parents.parent_idx == ep);
  OMEGA_H_CHECK(edge_parents.codes == ec);
  OMEGA_H_CHECK(face_parents.parent_idx == fp);
  OMEGA_H_CHECK(face_parents.codes == fc);
  OMEGA_H_CHECK(region_parents.parent_idx == rp);
  OMEGA_H_CHECK(region_parents.codes == rc);
}

static void check_3D_after(Omega_h::Mesh* m) {
  check_3D_leaves_after(m);
  check_3D_levels_after(m);
  check_3D_parents_after(m);
}

static void test_3D_arrays(Omega_h::Library* lib) {
  auto w = lib->world();
  auto f = OMEGA_H_HYPERCUBE;
  auto m = Omega_h::build_box(w, f, 2.0, 1.0, 1.0, 2, 1, 1);
  Int lengths[4] = {12, 20, 11, 2};
  check_before(&m, lengths);
  auto xfer_opts = Omega_h::TransferOpts();
  Omega_h::amr_refine(&m, Omega_h::Bytes({1, 0}), xfer_opts);
  check_3D_after(&m);
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  test_2D_arrays(&lib);
  test_3D_arrays(&lib);
}
