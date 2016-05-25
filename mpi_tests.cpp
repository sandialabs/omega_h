#include "internal.hpp"

static void test_one_rank(CommPtr comm) {
  CHECK(comm->size() == 1);
  { // make sure we can operate on zero-length data
  Dist dist;
  dist.set_comm(comm);
  dist.set_dest_ranks(Read<I32>({}));
  dist.set_dest_idxs(LOs({}), 0);
  dist.set_roots2items(LOs({0}));
  }
  {
  Dist dist;
  dist.set_comm(comm);
  dist.set_dest_ranks(Read<I32>({0,0,0,0}));
  dist.set_dest_idxs(LOs({3,2,1,0}), 4);
  Read<GO> a({0,1,2,3});
  auto b = dist.exch(a, 1);
  CHECK(b == Read<GO>({3,2,1,0}));
  }
}

static void test_two_ranks(CommPtr comm) {
  CHECK(comm->size() == 2);
  Dist dist;
  /* partition is {0,1,2}{3,4},
     global reversal to
                  {4,3,2}{1,0} */
  if (comm->rank() == 0) {
    dist.set_dest_ranks(Read<I32>({1,1,0}));
    dist.set_dest_idxs(LOs(       {1,0,2}), 3);
  } else {
    dist.set_dest_ranks(Read<I32>({0,0}));
    dist.set_dest_idxs(LOs(       {1,0}), 2);
  }
  Reals a;
  if (comm->rank() == 0) {
    a = Reals({0.,1.,2.});
  } else {
    a = Reals({3.,4.});
  }
  auto b = dist.exch(a, 1);
  if (comm->rank() == 0) {
    CHECK(b == Reals({4.,3.,2.}));
  } else {
    CHECK(b == Reals({1.,0.}));
  }
  auto c = dist.invert().exch(b, 1);
  CHECK(c == a);
}

static void test_all() {
  auto world = Comm::world();
  if (world->rank() == 0) {
    test_one_rank(Comm::self());
    auto one = world->split(world->rank(), 0);
    test_one_rank(one);
    test_one_rank(one->dup());
  }
  if (world->size() >= 2) {
    auto two = world->split(world->rank() / 2, world->rank() % 2);
    if (world->rank() / 2 == 0)
      test_two_ranks(two);
  }
}

int main(int argc, char** argv) {
  init(argc, argv);
  test_all();
  fini();
}
