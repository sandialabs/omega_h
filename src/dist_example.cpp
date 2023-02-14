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

    //Perform the exhcange operation, and verify it was executed correctly
    auto b = dist.exch(a, width);
    OMEGA_H_CHECK(b == Reals({2.0 * reverse_root + 1, 2.0 * reverse_root}));
}

static void test_reduce_exchange(CommPtr comm){
    OMEGA_H_CHECK(comm->size() == 3);
    Dist dist;
    dist.set_parent_comm(comm);
/*
 * Starting Partition: {1, 1}  {0, 0}  {0, 0}
 * End Partition:      {1, 1}  {1, 1}  {1, 1}
 * Reduced Exchange:   {3, 3}  {1, 1}  {1, 1}
*/
    if(comm->rank() == 0){
        dist.set_dest_ranks(Read<I32>({1, 2}));
        //Fan out using roots2items
        dist.set_roots2items(LOs({0, 2}));
    }
    else{
        dist.set_dest_ranks(Read<I32>({})); //Note: need to call on every rank
    }
    int width = 2;
    Reals a;

    if(comm->rank() == 0)
        a = Reals({1.0, 1.0});
    else
        a = Reals({});
    std::cout << comm->rank() << " Ready for Exchange" << std::endl; 
    auto b = dist.exch(a, width);
    if(comm->rank() == 0)
        OMEGA_H_CHECK(b == Reals({}));
    else{
	auto host_b = HostRead<Real>(b);
	std::cout << "Rank " << comm->rank() << ": Exchange Contents - " << host_b.get(0) 
		  << " " << host_b.get(1) << std::endl;	
        OMEGA_H_CHECK(b == Reals({1, 1}));
    }
}

int main(int argc, char** argv) {
    auto lib = Library(&argc, &argv);
    auto world = lib.world();
    assert(world->size() == 3);
    //test_three_ranks(world);
    test_reduce_exchange(world);
}
