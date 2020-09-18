#include <Omega_h_array_ops.hpp>
#include <Omega_h_bipart.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_inertia.hpp>
#include <Omega_h_owners.hpp>
#include <Omega_h_vtk.hpp>

#include <sstream>

using namespace Omega_h;

static void test_one_rank(CommPtr comm) {
  OMEGA_H_CHECK(comm->size() == 1);
  {  // make sure we can operate on zero-length data
    Dist dist;
    dist.set_parent_comm(comm);
    dist.set_dest_ranks(Read<I32>({}));
    dist.set_dest_idxs(LOs({}), 0);
    dist.set_roots2items(LOs({0}));
  }
  {
    Dist dist;
    dist.set_parent_comm(comm);
    dist.set_dest_ranks(Read<I32>({0, 0, 0, 0}));
    dist.set_dest_idxs(LOs({2, 3, 1, 0}), 4);
    //dist.set_dest_idxs(LOs({3, 2, 1, 0}), 4);
    Read<GO> a({0, 1, 2, 3});
    {
      auto b = dist.exch(a, 1);
      OMEGA_H_CHECK(b == Read<GO>({3, 2, 0, 1}));
      //OMEGA_H_CHECK(b == Read<GO>({3, 2, 1, 0}));
    }
  }
  {
    Dist copies2owners;
    copies2owners.set_parent_comm(comm);
    copies2owners.set_dest_ranks(Read<I32>({0, 0, 0, 0}));
    copies2owners.set_dest_idxs(LOs({0, 1, 2, 3}), 4);

    { // test variable sized array where actors have all size 1
      Read<GO> a({0, 1, 2, 3});
      auto const vdist = create_dist_for_variable_sized(copies2owners, {0, 1, 2, 3, 4});
      auto b = vdist.exch(a, 1);
      OMEGA_H_CHECK(b == Read<GO>({0, 1, 2, 3}));
    }
    { // test variable sized array where actors have size 1 but actors 2 and 3 which have size 2
      Read<GO> a({0, 1, 2, 3, 4, 5});
      auto const vdist = create_dist_for_variable_sized(copies2owners, {0, 1, 3, 4, 6});
      auto b = vdist.exch(a, 1);
      OMEGA_H_CHECK(b == Read<GO>({0, 1, 2, 3, 4, 5}));
    }
  }
}

static void test_two_ranks_dist(CommPtr comm) {
  OMEGA_H_CHECK(comm->size() == 2);
  Dist dist;
  dist.set_parent_comm(comm);
  /* partition is {0,1,2}{3,4},
     global reversal to
                  {4,3,2}{1,0} */
  if (comm->rank() == 0) {
    dist.set_dest_ranks(Read<I32>({1, 1, 0}));
    dist.set_dest_idxs(LOs({1, 0, 2}), 3);
  } else {
    dist.set_dest_ranks(Read<I32>({0, 0}));
    dist.set_dest_idxs(LOs({1, 0}), 2);
  }
  Reals a;
  if (comm->rank() == 0) {
    a = Reals({0., 1., 2.});
  } else {
    a = Reals({3., 4.});
  }
  auto b = dist.exch(a, 1);
  if (comm->rank() == 0) {
    OMEGA_H_CHECK(b == Reals({4., 3., 2.}));
  } else {
    OMEGA_H_CHECK(b == Reals({1., 0.}));
  }
  auto c = dist.invert().exch(b, 1);
  OMEGA_H_CHECK(c == a);
}

static void test_two_ranks_dist_for_two_variable_sized_actors(CommPtr comm) {
  OMEGA_H_CHECK(comm->size() == 2);
  Dist copies2owners;
  copies2owners.set_parent_comm(comm);
  /* partition is {0,1,2}{3,4},
     global reversal to
                  {3,4,2}{0,1} */
  if (comm->rank() == 0) {
    copies2owners.set_dest_ranks(Read<I32>({1, 1, 0}));
    copies2owners.set_dest_idxs(LOs({0, 1, 2}), 3);
  } else {
    copies2owners.set_dest_ranks(Read<I32>({0, 0}));
    copies2owners.set_dest_idxs(LOs({0, 1}), 2);
  }
  {
    LOs copies2data;
    if (comm->rank() == 0) {
      copies2data = {0, 1, 4, 6};
    } else {
      copies2data = {0, 1, 4};
    }
    auto dist = create_dist_for_variable_sized(copies2owners, copies2data);
    Reals a;
    if (comm->rank() == 0) {
      a = Reals({1, 3, 5, 7, 9, 11});
    } else {
      a = Reals({13, 15, 17, 19});
    }
    auto b = dist.exch(a, 1);
    if (comm->rank() == 0) {
      OMEGA_H_CHECK(b == Reals({13, 15, 17, 19, 9, 11}));
    } else {
      OMEGA_H_CHECK(b == Reals({1, 3, 5, 7}));
    }
    auto c = dist.invert().exch(b, 1);
    OMEGA_H_CHECK(c == a);
  }
  { // size of the first actor is 0
    LOs copies2data;
    if (comm->rank() == 0) {
      copies2data = {0, 0, 3, 5};
    } else {
      copies2data = {0, 0, 3};
    }
    auto dist = create_dist_for_variable_sized(copies2owners, copies2data);
    Reals a;
    if (comm->rank() == 0) {
      a = Reals({3, 5, 7, 9, 11});
    } else {
      a = Reals({15, 17, 19});
    }
    auto b = dist.exch(a, 1);
    if (comm->rank() == 0) {
      OMEGA_H_CHECK(b == Reals({15, 17, 19, 9, 11}));
    } else {
      OMEGA_H_CHECK(b == Reals({3, 5, 7}));
    }
    auto c = dist.invert().exch(b, 1);
    OMEGA_H_CHECK(c == a);
  }
  { // size of the second actor is 0
    LOs copies2data;
    if (comm->rank() == 0) {
      copies2data = {0, 1, 1, 3};
    } else {
      copies2data = {0, 1, 1};
    }
    auto dist = create_dist_for_variable_sized(copies2owners, copies2data);
    Reals a;
    if (comm->rank() == 0) {
      a = {1, 9, 11};
    } else {
      a = {13};
    }
    auto b = dist.exch(a, 1);
    if (comm->rank() == 0) {
      OMEGA_H_CHECK(b == Reals({13, 9, 11}));
    } else {
      OMEGA_H_CHECK(b == Reals({1}));
    }
    auto c = dist.invert().exch(b, 1);
    OMEGA_H_CHECK(c == a);
  }
}

static void test_two_rank_for_four_variable_sized_actors(CommPtr comm) {
  OMEGA_H_CHECK(comm->size() == 2);
  Dist copies2owners;
  copies2owners.set_parent_comm(comm);
  /* partition is {0,1,2,3}{4,5,6,7},
     global reversal to
                  {7,4,5,6}{1,2,3,0} */
  if (comm->rank() == 0) {
    copies2owners.set_dest_ranks(Read<I32>({1, 1, 1, 1}));
    copies2owners.set_dest_idxs(LOs({3, 0, 1, 2}), 4);
  } else {
    copies2owners.set_dest_ranks(Read<I32>({0, 0, 0, 0}));
    copies2owners.set_dest_idxs(LOs({1, 2, 3, 0}), 4);
  }
  LOs copies2data;
  if (comm->rank() == 0) {
    copies2data = LOs({0, 1, 4, 6, 7});
  } else {
    copies2data = LOs({0, 3, 5, 6, 7});
  }
  auto dist = create_dist_for_variable_sized(copies2owners, copies2data);
  Reals a;
  if (comm->rank() == 0) {
    a = Reals({3, -1, -3, -5, 5, 7, -7});
  } else {
    a = Reals({9, 11, 13, -9, -11, 15, -13});
  }
  auto b = dist.exch(a, 1);
  if (comm->rank() == 0) {
    OMEGA_H_CHECK(b == Reals({-13, 9, 11, 13, -9, -11, 15}));
  } else {
    OMEGA_H_CHECK(b == Reals({-1, -3, -5, 5, 7, -7, 3}));
  }
  auto c = dist.invert().exch(b, 1);
  OMEGA_H_CHECK(c == a);
}

static void test_two_ranks_eq_owners(CommPtr comm) {
  /* {0,1,2} and {2,3,4}.
     sizes are the same, so ownership of (2)
     will go to smallest rank, rank 0. */
  if (comm->rank() == 0) {
    auto owners = owners_from_globals(comm, Read<GO>({0, 1, 2}), Read<I32>());
    OMEGA_H_CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    OMEGA_H_CHECK(owners.idxs == Read<I32>({0, 1, 2}));
  } else {
    auto owners = owners_from_globals(comm, Read<GO>({2, 3, 4}), Read<I32>());
    OMEGA_H_CHECK(owners.ranks == Read<I32>({0, 1, 1}));
    OMEGA_H_CHECK(owners.idxs == Read<I32>({2, 1, 2}));
  }
}

static void test_two_ranks_uneq_owners(CommPtr comm) {
  /* {0,1,2} and {2,3}.
     now rank 1 has fewer items, so it
     will be deemed the owner */
  if (comm->rank() == 0) {
    auto owners = owners_from_globals(comm, Read<GO>({0, 1, 2}), Read<I32>());
    OMEGA_H_CHECK(owners.ranks == Read<I32>({0, 0, 1}));
    OMEGA_H_CHECK(owners.idxs == Read<I32>({0, 1, 0}));
  } else {
    auto owners = owners_from_globals(comm, Read<GO>({2, 3}), Read<I32>());
    OMEGA_H_CHECK(owners.ranks == Read<I32>({1, 1}));
    OMEGA_H_CHECK(owners.idxs == Read<I32>({0, 1}));
  }
}

static void test_two_ranks_forced_owners(CommPtr comm) {
  /* {0,1,2} and {2,3}.
     rank 1 still has fewer items, but we
     will force 0 to be the owner by specifying
     owner ranks for the copies */
  if (comm->rank() == 0) {
    auto owners =
        owners_from_globals(comm, Read<GO>({0, 1, 2}), Read<I32>({0, 0, 0}));
    OMEGA_H_CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    OMEGA_H_CHECK(owners.idxs == Read<I32>({0, 1, 2}));
  } else {
    auto owners =
        owners_from_globals(comm, Read<GO>({2, 3}), Read<I32>({0, 1}));
    OMEGA_H_CHECK(owners.ranks == Read<I32>({0, 1}));
    OMEGA_H_CHECK(owners.idxs == Read<I32>({2, 1}));
  }
}

static void test_two_ranks_bipart(CommPtr comm) {
  if (comm->rank() == 0) {
    auto dist = bi_partition(comm, Read<I8>{0, 1, 0, 1});
    auto out = dist.exch(LOs({0, 1, 2, 3}), 1);
    OMEGA_H_CHECK(out == LOs({0, 2, 4}));
  } else {
    auto dist = bi_partition(comm, Read<I8>{0, 1, 1});
    auto out = dist.exch(LOs({4, 5, 6}), 1);
    OMEGA_H_CHECK(out == LOs({1, 3, 5, 6}));
  }
}

static void test_two_ranks_exch_sum(CommPtr comm) {
  Dist dist;
  dist.set_parent_comm(comm);
  /* three copies on each side, two are shared and
     owned by rank 0 */
  if (comm->rank() == 0) {
    dist.set_dest_ranks(Read<I32>({0, 0, 0}));
    dist.set_dest_idxs(LOs({0, 1, 2}), 3);
  } else {
    dist.set_dest_ranks(Read<I32>({0, 0, 1}));
    dist.set_dest_idxs(LOs({1, 2, 0}), 1);
  }
  auto recvd = dist.exch_reduce(LOs({1, 1, 1}), 1, OMEGA_H_SUM);
  if (comm->rank() == 0) {
    OMEGA_H_CHECK(recvd == LOs({1, 2, 2}));
  } else {
    OMEGA_H_CHECK(recvd == LOs({1}));
  }
}

static void test_two_ranks_owners(CommPtr comm) {
  test_two_ranks_eq_owners(comm);
  test_two_ranks_uneq_owners(comm);
  test_two_ranks_forced_owners(comm);
}

static void test_resolve_derived(CommPtr comm) {
  auto verts2globs = Read<GO>();
  if (comm->rank() == 0) {
    verts2globs = Read<GO>({0, 1, 3});
  } else {
    verts2globs = Read<GO>({1, 2, 3});
  }
  auto ev2v = LOs({0, 1, 1, 2, 2, 0});
  auto owners = Remotes();
  resolve_derived_copies(comm, verts2globs, 2, &ev2v, &owners);
  /* in both cases, the last edge is flipped since it goes
     from high global number to low. */
  OMEGA_H_CHECK(ev2v == LOs({0, 1, 1, 2, 0, 2}));
  if (comm->rank() == 0) {
    OMEGA_H_CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    OMEGA_H_CHECK(owners.idxs == Read<LO>({0, 1, 2}));
  } else {
    OMEGA_H_CHECK(owners.ranks == Read<I32>({1, 1, 0}));
    OMEGA_H_CHECK(owners.idxs == Read<LO>({0, 1, 1}));
  }
}

static void test_construct(Library* lib, CommPtr comm) {
  auto verts2globs = Read<GO>();
  if (comm->rank() == 0) {
    verts2globs = Read<GO>({0, 1, 3});
  } else {
    verts2globs = Read<GO>({1, 2, 3});
  }
  auto tv2v = LOs({0, 1, 2});
  Mesh mesh(lib);
  build_from_elems2verts(&mesh, comm, OMEGA_H_SIMPLEX, 2, tv2v, verts2globs);
  auto ev2v = mesh.ask_verts_of(EDGE);
  auto owners = mesh.ask_owners(EDGE);
  OMEGA_H_CHECK(ev2v == LOs({0, 1, 0, 2, 1, 2}));
  if (comm->rank() == 0) {
    OMEGA_H_CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    OMEGA_H_CHECK(owners.idxs == Read<LO>({0, 1, 2}));
  } else {
    OMEGA_H_CHECK(owners.ranks == Read<I32>({1, 0, 1}));
    OMEGA_H_CHECK(owners.idxs == Read<LO>({0, 2, 2}));
  }
}
/*
static void test_read_vtu(Library* lib, CommPtr comm) {
  auto mesh0 = build_box(comm, OMEGA_H_SIMPLEX, 1., 1., 0., 1, 1, 0);
  std::stringstream stream;
  vtk::write_vtu(
      stream, &mesh0, mesh0.dim(), vtk::get_all_vtk_tags(&mesh0, mesh0.dim()));
  Mesh mesh1(lib);
  vtk::read_vtu(stream, comm, &mesh1);
  auto opts = MeshCompareOpts::init(&mesh0, VarCompareOpts::zero_tolerance());
  OMEGA_H_CHECK(
      OMEGA_H_SAME == compare_meshes(&mesh0, &mesh1, opts, true, false));
}
*/
static void test_binary_io(Library* lib, CommPtr comm) {
  auto mesh0 = build_box(comm, OMEGA_H_SIMPLEX, 1., 1., 0., 4, 4, 0);
  mesh0.set_parting(OMEGA_H_ELEM_BASED);
  binary::write("mpi_test_elem_based.osh", &mesh0);
  Mesh mesh1(lib);
  binary::read("mpi_test_elem_based.osh", comm, &mesh1);
  auto opts = MeshCompareOpts::init(&mesh0, VarCompareOpts::zero_tolerance());
  OMEGA_H_CHECK(
      OMEGA_H_SAME == compare_meshes(&mesh0, &mesh1, opts, true, true));
  mesh0.set_parting(OMEGA_H_GHOSTED);
  binary::write("mpi_test_ghosted.osh", &mesh0);
  Mesh mesh2(lib);
  binary::read("mpi_test_ghosted.osh", comm, &mesh2);
  OMEGA_H_CHECK(
      OMEGA_H_SAME == compare_meshes(&mesh0, &mesh2, opts, true, true));
}

static void test_two_ranks(Library* lib, CommPtr comm) {
  test_two_ranks_dist(comm);
  test_two_ranks_dist_for_two_variable_sized_actors(comm);
  test_two_rank_for_four_variable_sized_actors(comm);
  test_two_ranks_owners(comm);
  test_two_ranks_bipart(comm);
  test_two_ranks_exch_sum(comm);
  test_resolve_derived(comm);
  test_construct(lib, comm);
  //test_read_vtu(lib, comm);
  test_binary_io(lib, comm);
}

void test_rib(CommPtr comm) {
  auto rank = comm->rank();
  auto size = comm->size();
  LO n = 5;
  Write<Real> w_coords(n * 3);
  auto set_coords = OMEGA_H_LAMBDA(LO i) {
    set_vector(w_coords, i, vector_3(i * size + rank, 0, 0));
  };
  parallel_for(n, set_coords);
  Reals coords(w_coords);
  Reals masses(n, 1);
  auto owners = Remotes(Read<I32>(n, rank), LOs(n, 0, 1));
  auto hints = inertia::Rib();
  inertia::recursively_bisect(comm, 1.1, &coords, &masses, &owners, &hints);
  I32 size2 = 1;
  for (auto axis : hints.axes) {
    OMEGA_H_CHECK(are_close(axis, vector_3(1, 0, 0)));
    size2 *= 2;
  }
  OMEGA_H_CHECK(size2 == size);
  auto check_coords = OMEGA_H_LAMBDA(LO i) {
    auto v = get_vector<3>(coords, i);
    OMEGA_H_CHECK(rank * n <= v[0]);
    OMEGA_H_CHECK(v[0] < (rank + 1) * n);
    OMEGA_H_CHECK(v[1] == 0 && v[2] == 0);
  };
  parallel_for(n, check_coords);
  OMEGA_H_CHECK(masses == Reals(n, 1));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  if (world->rank() == 0) {
    test_one_rank(lib.self());
  }
  auto one = world->split(world->rank(), 0);
  if (world->rank() == 0) {
    test_one_rank(one);
    test_one_rank(one->dup());
  }
  if (world->size() >= 2) {
    auto two = world->split(world->rank() / 2, world->rank() % 2);
    if (world->rank() / 2 == 0) {
      test_two_ranks(&lib, two);
    }
  }
  world->barrier();
  test_rib(world);
}
