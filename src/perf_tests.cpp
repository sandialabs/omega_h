#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <random>

#include "access.hpp"
#include "adjacency.hpp"
#include "array.hpp"
#include "eigen.hpp"
#include "internal.hpp"
#include "loop.hpp"
#include "metric.hpp"
#include "sort.hpp"
#include "space.hpp"
#include "timer.hpp"

using namespace Omega_h;

static Int const nelems = 1000 * 1000;

/* Intel compiler version 16 doesn't seem to have
 * std::uniform_real_distribution and std::uniform_int_distribution */
#ifndef __INTEL_COMPILER

static Reals random_reals(Int n, Real from, Real to) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(from, to);
  HostWrite<Real> ha(n);
  for (Int i = 0; i < n; ++i) ha[i] = dis(gen);
  return ha.write();
}

/* we'll generate metric tensors which have random
   rotation matrices (choose two random Euler angles.
   this is not uniform over the sphere, but thats not important),
   and desired lengths (1,1,a^2), where (a) is the desired
   amount of anisotropy (a=1000 implies 1000:1 ratio). */
static double const anisotropy = 1000;

static Reals random_metrics() {
  Reals alphas = random_reals(nelems, 0, PI / 2);
  Reals betas = random_reals(nelems, 0, PI / 2);
  Write<Real> write_metrics(nelems * 6);
  auto f0 = LAMBDA(Int i) {
    auto r = rotate(alphas[i], vector_3(0, 0, 1)) *
             rotate(betas[i], vector_3(0, 1, 0));
    auto m = compose_metric(r, vector_3(1., 1., 1.0 / anisotropy));
    set_symm(write_metrics, i, m);
  };
  parallel_for(nelems, f0);
  return write_metrics;
}

static void test_metric_decompose(Reals metrics) {
  /* now, decompose the metrics and get the largest
     eigenvalue of each */
  Write<Real> write_eigenvs(nelems);
  auto f1 = LAMBDA(Int i) {
    auto m = get_symm<3>(metrics, i);
    Matrix<3, 3> r;
    Vector<3> l;
    decompose_eigen(m, r, l);
    auto eigenv = max2(max2(l[0], l[1]), l[2]);
    write_eigenvs[i] = eigenv;
  };
  Now t0 = now();
  Int niters = 3;
  for (Int i = 0; i < niters; ++i) parallel_for(nelems, f1);
  Now t1 = now();
  std::cout << "eigendecomposition of " << nelems << " metric tensors "
            << niters << " times takes " << (t1 - t0) << " seconds\n";
  CHECK(are_close(Reals(write_eigenvs), Reals(nelems, square(anisotropy))));
}

static void test_metric_invert(Reals metrics) {
  /* now, decompose the metrics and get the largest
     eigenvalue of each */
  Write<Real> write_vals(nelems);
  auto f1 = LAMBDA(Int i) {
    auto m = get_symm<3>(metrics, i);
    auto inv = invert(m);
    write_vals[i] = max_norm(inv);
  };
  Now t0 = now();
  Int niters = 30;
  for (Int i = 0; i < niters; ++i) parallel_for(nelems, f1);
  Now t1 = now();
  std::cout << "inversion of " << nelems << " metric tensors " << niters
            << " times takes " << (t1 - t0) << " seconds\n";
}

static void test_metric_math() {
  Reals metrics = random_metrics();
  test_metric_decompose(metrics);
  test_metric_invert(metrics);
}

unsigned uniform(
    unsigned m); /* Returns a random integer 0 <= uniform(m) <= m-1 */

/* Fisher-Yates shuffle, a.k.a Knuth shuffle */
static Read<Int> random_perm(Int n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Int> dis;
  /* start with identity permutation */
  HostWrite<Int> permutation(n, 0, 1);
  for (Int i = 0; i + 1 < n; i++) {
    Int j = dis(gen) % (n - i); /* A random integer such that 0 ≤ j < n-i*/
    std::swap(permutation[i], permutation[i + j]);
  }
  return Read<Int>(permutation.write());
}

static void test_repro_sum() {
  Reals inputs = random_reals(nelems, 0, 1e100);
  Real rs = 0, s = 0;
  Int niters = 100;
  {
    Now t0 = now();
    for (Int i = 0; i < niters; ++i) rs = repro_sum(inputs);
    Now t1 = now();
    std::cout << "reproducibly adding " << nelems << " reals " << niters
              << " times "
              << "takes " << (t1 - t0) << " seconds\n";
  }
  {
    Now t0 = now();
    for (Int i = 0; i < niters; ++i) s = sum(inputs);
    Now t1 = now();
    std::cout << "adding " << nelems << " reals " << niters << " times "
              << "takes " << (t1 - t0) << " seconds\n";
  }
  CHECK(are_close(s, rs));
  Read<Int> p = random_perm(nelems);
  Write<Real> write_shuffled(nelems);
  auto f = LAMBDA(Int i) { write_shuffled[i] = inputs[p[i]]; };
  parallel_for(nelems, f);
  Reals shuffled(write_shuffled);
  Real rs2 = repro_sum(shuffled);
  Real s2 = sum(shuffled);
  CHECK(are_close(s2, rs2));
  CHECK(rs == rs2); /* bitwise reproducibility ! */
  if (s == s2) std::cerr << "warning: the naive sum gave the same answer\n";
}

template <typename T>
static Read<T> random_ints(Int n, T from, T to) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<T> dis(from, to);
  HostWrite<T> ha(n);
  for (Int i = 0; i < n; ++i) ha[i] = dis(gen);
  return ha.write();
}

static void test_sort_n(Int width) {
  LOs a = random_ints<LO>(nelems * width, 0, nelems);
  LOs perm;
  Int niters = 5;
  Now t0 = now();
  for (Int i = 0; i < niters; ++i) perm = sort_by_keys(a, width);
  Now t1 = now();
  std::cout << "sorting " << nelems << " sets of " << width << " integers "
            << niters << " times takes " << (t1 - t0) << " seconds\n";
}

static void test_sort() {
  test_sort_n(1);
  test_sort_n(2);
  test_sort_n(3);
}

#endif

static void test_invert_adj(LOs tets2verts, LO nverts) {
  LO ntets = tets2verts.size() / 4;
  Read<GO> tet_globals(ntets, 0, 1);
  Adj inv;
  Int niters = 5;
  Now t0 = now();
  for (Int i = 0; i < niters; ++i)
    inv = invert_adj(Adj(tets2verts), 4, nverts, tet_globals);
  Now t1 = now();
  std::cout << "inverting " << ntets << " tets -> verts " << niters
            << " times takes " << (t1 - t0) << " seconds\n";
}

static void test_reflect_down(LOs tets2verts, LOs tris2verts, LO nverts) {
  LO ntets = tets2verts.size() / 4;
  LO ntris = tris2verts.size() / 3;
  Int niters = 2;
  Adj verts2tris =
      invert_adj(Adj(tris2verts), 3, nverts, Read<GO>(ntris, 0, 1));
  {
    Now t0 = now();
    for (Int i = 0; i < niters; ++i)
      reflect_down(tets2verts, tris2verts, verts2tris, 3, 2);
    Now t1 = now();
    std::cout << "reflect_down " << ntets << " tets -> tris "
              << "by only upward " << niters << " times takes " << (t1 - t0)
              << " seconds\n";
  }
}

static void test_adjs(Library const& lib) {
  Mesh mesh;
  {
    Now t0 = now();
    auto nx = 42;
    build_box(&mesh, lib, 1, 1, 1, nx, nx, nx);
    Now t1 = now();
    std::cout << "building a " << nx << 'x' << nx << 'x' << nx << " box took "
              << (t1 - t0) << " seconds\n";
  }
  {
    Now t0 = now();
    mesh.reorder();
    Now t1 = now();
    std::cout << "reordering a " << mesh.nelems() << " tet mesh took "
              << (t1 - t0) << " seconds\n";
  }
  LOs tets2verts;
  LOs tris2verts;
  {
    Now t0 = now();
    tets2verts = mesh.ask_verts_of(TET);
    tris2verts = mesh.ask_verts_of(TRI);
    Now t1 = now();
    std::cout << "asking tet->vert and tri->vert of a " << mesh.nelems()
              << " tet mesh took " << (t1 - t0) << " seconds\n";
  }
  auto nverts = mesh.nverts();
  test_invert_adj(tets2verts, nverts);
  test_reflect_down(tets2verts, tris2verts, nverts);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
#ifndef __INTEL_COMPILER
  test_metric_math();
  test_repro_sum();
  test_sort();
#endif
  test_adjs(lib);
}
