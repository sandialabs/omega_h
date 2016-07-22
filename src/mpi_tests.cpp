#include "all.hpp"

using namespace osh;

static void test_one_rank(CommPtr comm) {
  CHECK(comm->size() == 1);
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
    dist.set_dest_idxs(LOs({3, 2, 1, 0}), 4);
    Read<GO> a({0, 1, 2, 3});
    auto b = dist.exch(a, 1);
    CHECK(b == Read<GO>({3, 2, 1, 0}));
  }
}

static void test_two_ranks_dist(CommPtr comm) {
  CHECK(comm->size() == 2);
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
    CHECK(b == Reals({4., 3., 2.}));
  } else {
    CHECK(b == Reals({1., 0.}));
  }
  auto c = dist.invert().exch(b, 1);
  CHECK(c == a);
}

static void test_two_ranks_eq_owners(CommPtr comm) {
  /* {0,1,2} and {2,3,4}.
     sizes are the same, so ownership of (2)
     will go to smallest rank, rank 0. */
  if (comm->rank() == 0) {
    auto owners = owners_from_globals(comm, Read<GO>({0, 1, 2}), Read<I32>());
    CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    CHECK(owners.idxs == Read<I32>({0, 1, 2}));
  } else {
    auto owners = owners_from_globals(comm, Read<GO>({2, 3, 4}), Read<I32>());
    CHECK(owners.ranks == Read<I32>({0, 1, 1}));
    CHECK(owners.idxs == Read<I32>({2, 1, 2}));
  }
}

static void test_two_ranks_uneq_owners(CommPtr comm) {
  /* {0,1,2} and {2,3}.
     now rank 1 has fewer items, so it
     will be deemed the owner */
  if (comm->rank() == 0) {
    auto owners = owners_from_globals(comm, Read<GO>({0, 1, 2}), Read<I32>());
    CHECK(owners.ranks == Read<I32>({0, 0, 1}));
    CHECK(owners.idxs == Read<I32>({0, 1, 0}));
  } else {
    auto owners = owners_from_globals(comm, Read<GO>({2, 3}), Read<I32>());
    CHECK(owners.ranks == Read<I32>({1, 1}));
    CHECK(owners.idxs == Read<I32>({0, 1}));
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
    CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    CHECK(owners.idxs == Read<I32>({0, 1, 2}));
  } else {
    auto owners =
        owners_from_globals(comm, Read<GO>({2, 3}), Read<I32>({0, 1}));
    CHECK(owners.ranks == Read<I32>({0, 1}));
    CHECK(owners.idxs == Read<I32>({2, 1}));
  }
}

static void test_two_ranks_bipart(CommPtr comm) {
  if (comm->rank() == 0) {
    auto dist = bi_partition(comm, Read<I8>{0, 1, 0, 1});
    auto out = dist.exch(LOs({0, 1, 2, 3}), 1);
    CHECK(out == LOs({0, 2, 4}));
  } else {
    auto dist = bi_partition(comm, Read<I8>{0, 1, 1});
    auto out = dist.exch(LOs({4, 5, 6}), 1);
    CHECK(out == LOs({1, 3, 5, 6}));
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
  auto recvd = dist.exch_reduce(LOs({1, 1, 1}), 1, OSH_SUM);
  if (comm->rank() == 0) {
    CHECK(recvd == LOs({1, 2, 2}));
  } else {
    CHECK(recvd == LOs({1}));
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
  CHECK(ev2v == LOs({0, 1, 1, 2, 0, 2}));
  if (comm->rank() == 0) {
    CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    CHECK(owners.idxs == Read<LO>({0, 1, 2}));
  } else {
    CHECK(owners.ranks == Read<I32>({1, 1, 0}));
    CHECK(owners.idxs == Read<LO>({0, 1, 1}));
  }
}

static void test_construct(CommPtr comm) {
  auto verts2globs = Read<GO>();
  if (comm->rank() == 0) {
    verts2globs = Read<GO>({0, 1, 3});
  } else {
    verts2globs = Read<GO>({1, 2, 3});
  }
  auto tv2v = LOs({0, 1, 2});
  Mesh mesh;
  build_from_elems2verts(&mesh, comm, 2, tv2v, verts2globs);
  auto ev2v = mesh.ask_verts_of(EDGE);
  auto owners = mesh.ask_owners(EDGE);
  CHECK(ev2v == LOs({0, 1, 0, 2, 1, 2}));
  if (comm->rank() == 0) {
    CHECK(owners.ranks == Read<I32>({0, 0, 0}));
    CHECK(owners.idxs == Read<LO>({0, 1, 2}));
  } else {
    CHECK(owners.ranks == Read<I32>({1, 0, 1}));
    CHECK(owners.idxs == Read<LO>({0, 2, 2}));
  }
}

static void test_read_vtu(Library const& lib, CommPtr comm) {
  Mesh mesh0;
  if (comm->rank() == 0) {
    build_box(&mesh0, lib, 1, 1, 0, 1, 1, 0);
  }
  mesh0.set_comm(comm);
  mesh0.balance();
  std::stringstream stream;
  vtk::write_vtu(stream, &mesh0, mesh0.dim());
  Mesh mesh1;
  vtk::read_vtu(stream, comm, &mesh1);
  CHECK(OSH_SAME == compare_meshes(&mesh0, &mesh1, 0.0, 0.0, true, false));
}

static void test_two_ranks(Library const& lib, CommPtr comm) {
  test_two_ranks_dist(comm);
  test_two_ranks_owners(comm);
  test_two_ranks_bipart(comm);
  test_two_ranks_exch_sum(comm);
  test_resolve_derived(comm);
  test_construct(comm);
  test_read_vtu(lib, comm);
}

static void test_rib(CommPtr comm) {
  auto rank = comm->rank();
  auto size = comm->size();
  LO n = 5;
  Write<Real> w_coords(n * 3);
  auto set_coords = LAMBDA(LO i) {
    set_vector(w_coords, i, vector_3(i * size + rank, 0, 0));
  };
  parallel_for(n, set_coords);
  Reals coords(w_coords);
  Reals masses(n, 1);
  auto owners = Remotes(Read<I32>(n, rank), LOs(n, 0, 1));
  auto out = inertia::recursively_bisect(comm, coords, masses, owners, 0.01,
                                         inertia::Rib());
  I32 size2 = 1;
  for (auto axis : out.axes) {
    CHECK(are_close(axis, vector_3(1, 0, 0)));
    size2 *= 2;
  }
  CHECK(size2 == size);
  auto check_coords = LAMBDA(LO i) {
    auto v = get_vector<3>(coords, i);
    CHECK(rank * n <= v[0]);
    CHECK(v[0] < (rank + 1) * n);
    CHECK(v[1] == 0 && v[2] == 0);
  };
  parallel_for(n, check_coords);
  CHECK(masses == Reals(n, 1));
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
      test_two_ranks(lib, two);
    }
  }
  test_rib(world);
}
