#include "internal.hpp"

static void test_one_rank(CommPtr comm) {
  CHECK(comm->size() == 1);
  Dist dist;
  dist.set_comm(comm);
  dist.set_dest_ranks(Read<I32>({0,0,0,0}));
  dist.set_dest_idxs(Read<I32>({3,2,1,0}), 4);
  Read<GO> a({0,1,2,3});
  auto b = dist.exch(a, 1);
  CHECK(b == Read<GO>({3,2,1,0}));
}

static void test_all() {
  auto world = Comm::world();
  if (world->rank() == 0) {
    test_one_rank(Comm::self());
    auto one = world->split(world->rank(), 0);
    test_one_rank(one);
    test_one_rank(one->dup());
  }
  if (world->size() == 1)
    test_one_rank(world);
}

int main(int argc, char** argv) {
  init(argc, argv);
  test_all();
  fini();
}
