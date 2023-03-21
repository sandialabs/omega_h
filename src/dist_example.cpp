#include <Omega_h_array_ops.hpp>
#include <Omega_h_dist.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_library.hpp>

using namespace Omega_h;

static void test_three_ranks(CommPtr comm){
  OMEGA_H_CHECK(comm->size() == 3);
  Dist dist;
  dist.set_parent_comm(comm);
/*
 * Starting Partition: {0, 1}  {2, 3}  {4, 5}
 * Reversal to       : {5, 4}  {3, 2}  {1, 0}
*/
  //Configure the directed graph by finding our destination node(reverse root)
  int reverse_root = 2 - comm->rank();
  dist.set_dest_ranks(Read<I32>({reverse_root, reverse_root}));
/*
 * Forward: Rank 0    Rank 1    Rank 2
 *            │ ┏━━━━━━╾╫╼━━━━━━━┛                 
 *            ╰─╂───────╫─────────╮
 *            ┏━┛       ║         │
 * Reverse: Rank 0    Rank 1    Rank 2
*/ 
  //Set the destination indexes of the reverse root
  dist.set_dest_idxs(LOs({1, 0}), 2);

  //Generate the Packets for the forward root(2 Packets of width 1)
  int width = 1;
  Reals a({2.0 * comm->rank(), 2.0 * comm->rank() + 1});

  //Perform the exchange operation, and verify it was executed correctly
  auto b = dist.exch(a, width);
  OMEGA_H_CHECK(b == Reals({2.0 * reverse_root + 1, 2.0 * reverse_root}));
}

static void test_packet_fanout(CommPtr comm){
  OMEGA_H_CHECK(comm->size() == 3);
  Dist dist;
  dist.set_parent_comm(comm);
/*
 * Starting Partition: {1, 2}  {0, 0}  {0, 0}
 * End Partition:      {}      {1, 1}  {2, 2}
*/

  //Configure the directed graph. Note when explicitly setting parameters such as roots2items, 
  //undefined behavior may occur when not set for every rank. Also note that roots2items is an
  //offset array.
  if(comm->rank() == 0){
    dist.set_dest_ranks(Read<I32>({1, 1, 2, 2}));
    dist.set_roots2items(LOs({0, 2, 4}));
  }
  else{
    dist.set_dest_ranks(Read<I32>({})); 
    //Offset arrays of length 0 may cause undefined behavior, especially in inverted dists
    dist.set_roots2items(LOs({0}));
  }

  //Allocate Starting Partition
  Reals a;
  if(comm->rank() == 0)
    a = Reals({1.0, 2.0});
  else
    a = Reals({});
  
  //Perform Exchange Operation
  int width = 1;
  auto b = dist.exch(a, width);
  
  //Check our that our end partition is what we expected
  if(comm->rank() == 0)
    OMEGA_H_CHECK(b == Reals({}));
  else if(comm->rank() == 1){
    OMEGA_H_CHECK(b == Reals({1, 1}));
  }
  else{
    OMEGA_H_CHECK(b == Reals({2, 2}));
  }
}

static void test_rank_fanout(CommPtr comm){
  OMEGA_H_CHECK(comm->size() == 3);
  Dist dist;
  dist.set_parent_comm(comm);
/*
 * Starting Partition: {1, 2}  {0, 0}  {0, 0}
 * End Partition:      {}      {1, 2}  {1, 2}
*/

  //Configure the directed graph. Note when explicitly setting parameters such as roots2items, 
  //undefined behavior may occur when not set for every rank. Also note that roots2items is an
  //offset array.
  if(comm->rank() == 0){
    dist.set_dest_ranks(Read<I32>({1, 2, 1, 2}));
    dist.set_roots2items(LOs({0, 2, 4}));
  }
  else{
    dist.set_dest_ranks(Read<I32>({}));
    //Offset arrays of length 0 may cause undefined behavior, especially in inverted dists
    dist.set_roots2items(LOs({0}));
  }

  //Allocate Starting Partition
  Reals a;
  if(comm->rank() == 0)
    a = Reals({1.0, 2.0});
  else
    a = Reals({});

  //Perform Exchange Operation
  int width = 1;
  auto b = dist.exch(a, width);

  //Check End Partition is what we expected
  if(comm->rank() == 0)
    OMEGA_H_CHECK(b == Reals({}));
  else{
    OMEGA_H_CHECK(b == Reals({1, 2}));
  }
}

static void test_exchange_reduce(CommPtr comm){
  OMEGA_H_CHECK(comm->size() == 3);
  Dist dist;
  dist.set_parent_comm(comm);
/*
 * Starting Partition: {}   {1}  {1}
 * End Partition:      {2}  {}   {}
*/
  
  //Configure the directed graph. Note that dest_idxs must be set for all destination ranks, 
  //as well as the number of destination roots on this rank.
  if(comm->rank() == 0){
    dist.set_dest_ranks(Read<I32>({}));
    dist.set_dest_idxs(LOs({}), 1);
  }
  else{
    dist.set_dest_ranks(Read<I32>({0}));
    dist.set_dest_idxs(LOs({0}), 0);
  }

  //Allocate Starting Position
  Reals a;
  if(comm->rank() == 0)
    a = Reals({});
  else
    a = Reals({1});
 
  //Perform Exchange Operation
  int width = 1;
  auto b = dist.exch_reduce(a, width, OMEGA_H_SUM);	

  //Check End Partition is what we expected
  if(comm->rank() == 0)
    OMEGA_H_CHECK(b == Reals({2}));
  else
    OMEGA_H_CHECK(b == Reals({}));
}

static void test_fan_reduce(CommPtr comm){
  OMEGA_H_CHECK(comm->size() == 3);
  Dist dist;
  dist.set_parent_comm(comm);
/*
 * Starting Partition: {1}  {}   {}
 * End Partition:      {}   {1}  {1}
 * Inverted Partition: {2}  {}   {}
*/

  //Configure the directed graph. Note when explicitly setting parameters such as roots2items, 
  //undefined behavior may occur when not set for every rank. Also note that roots2items is an
  //offset array.
  if(comm->rank() == 0){
    dist.set_dest_ranks(Read<I32>({1, 2}));
    dist.set_roots2items(LOs({0, 2}));
    }
  else{
    dist.set_dest_ranks(Read<I32>({}));
    dist.set_roots2items(LOs({0}));
  }

  //Allocate Starting Partition
  Reals a;
  if(comm->rank() == 0)
    a = Reals({1});
  else
    a = Reals({});

  //Perform Exchange Operation		
  int width = 1;
  auto b = dist.exch(a, width);	

  //Check End Partition is what we expected
  if(comm->rank() == 0)
    OMEGA_H_CHECK(b == Reals({}));
  else
    OMEGA_H_CHECK(b == Reals({1}));

  //Invert Dist and Perform Exchange Reduce Operation
  //Note: The Dist Class stores 2 copies of roots2items, one forward and one reverse. It is critical
  //that roots2items is set to an array of non-zero length before the inversion. Failure to do so will 
  //result in undefined behavior
  Dist invDist = dist.invert();
  auto c = invDist.exch_reduce(b, width, OMEGA_H_SUM);

  //Check Inverted Partition is what we expect
  if(comm->rank() == 0)
    OMEGA_H_CHECK(c == Reals({2}));
  else
    OMEGA_H_CHECK(c == Reals({}));
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto world = lib.world();
  assert(world->size() == 3);
  test_three_ranks(world);
  test_packet_fanout(world);
  test_rank_fanout(world);
  test_exchange_reduce(world);
  test_fan_reduce(world);
}
