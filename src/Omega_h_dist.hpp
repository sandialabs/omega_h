#ifndef OMEGA_H_DIST_HPP
#define OMEGA_H_DIST_HPP

#include <Omega_h_comm.hpp>
#include <Omega_h_remotes.hpp>
#include <Omega_h_graph.hpp>

/* \file Omega_h_dist.hpp
  \brief Header file for Omega_h::Dist
 */

namespace Omega_h {

/*! \brief an MPI distributor object which encapsulates the idea
          of a communication pattern between lots of small
          actors on each MPI rank.
   
  \details
   Welcome to Dist, the magical parallel machine !

   This class implements a communication pattern
   that travels along edges of a partitioned graph.
   The graph is bipartite, with the two sets of nodes
   being the "source" set and the "destination" set.
   Both the source set and the destination set are
   partitioned among MPI ranks, such that each graph
   node resides on exactly one MPI rank.
   Arbitrary graph edges may be defined, i.e. each
   source may send to several destinations,
   and each destination may receive from several sources.
   There is a beautiful symmetry here that we take advantage
   of, namely that reversing the graph by swapping the
   source and destination sets produces a very similar structure.
   Dist is designed primarily to carry data from the
   source set to the destination set, but can also
   derive and provide much useful information about the
   parallel graph itself.

   The communication "packets" that Dist handles are
   tuples of N values of type T, where T is one of the
   basic types (integer, real).
   For efficiency, all packets in one MPI rank are
   packed into a single array.

   By performing several transformations on this array,
   Dist is able to perform arbitrary communication patterns.

   The array begins in a form that is sorted by
   source graph nodes on the local MPI rank.
   We call graph nodes on one rank "roots".
   "Forward roots" are the source graph nodes on the local MPI rank.
   So the array begins as one packet (tuple) per forward root.

   Then the array is expanded by duplicating packets until
   it has one packet per graph edge whose source node is on
   this MPI rank. At that stage, the graph edges are sorted
   by source node.
   We call graph edges sorted by graph node "items".

   Our MPI communication workhorse is the noble and tireless
   MPI_Neighbor_alltoallv call introduced in the MPI 3.0 standard.
   This function requires its input and output arrays to be
   sorted by MPI rank.
   As such, there is a permutation from "items" (graph edges
   sorted by source graph node) to "content" (graph edges sorted
   by MPI rank of destination graph node).

   Once the array is in "content" form, it is given to
   MPI_Neighbor_alltoallv which sends the contents around,
   and returns an array (also in "content" form) of received
   data.

   The received data follows the reverse path of the sent data.
   It is first "unpermuted" from "content" form to "items" form.

   At this point, we have a choice to make.
   Going from "items" to "roots" involves combining several
   packets together into one (reduction).
   Although the user may want a standard reduction (such as
   a sum), in many cases in Omega_h, more complicated operations
   need to be performed on this data.
   As such, the default API (exch()) will stop at destination
   items form and return the received data there.
   The exch_reduce() API will perform one of the standard
   pre-defined reductions, returning the array in "roots" form.

   Due to the symmetry, we use the terms "forward" and "reverse"
   to refer to roots, items, content, and communicators for the
   sent and received data, respectively.
 */

class Dist {
  CommPtr parent_comm_;
  LOs roots2items_[2];
  LOs items2content_[2];
  LOs msgs2content_[2];
  CommPtr comm_[2];

 public:
  Dist();
  Dist(Dist const& other);
  Dist& operator=(Dist const& other);
  Dist(CommPtr comm_in, Remotes fitems2rroots, LO nrroots);
  void set_parent_comm(CommPtr parent_comm_in);
  /* set the destination graph node MPI ranks of the "items"
     (graph edges sorted by local source node).
     this call defines the majority of the communication pattern,
     including neighbors and forwards items2content permutation.
     this is a collective call.
   */
  void set_dest_ranks(Read<I32> items2ranks_in);
  /* optionally set the destination root (local) indices of
     each forward item. one also needs to specify the number
     of destination roots on this MPI rank.
     this will specify the receiving side of the data
     structure, including items2content and roots2items.
     if not specified, then destination content, items,
     and roots will all be assumed to be the same.
   */
  void set_dest_idxs(LOs fitems2rroots, LO nrroots);
  /* optionally specify the expansion ("fan") from forward
     roots to items.
     if not specified, forward roots and items will be assumed
     to be the same. */
  void set_roots2items(LOs froots2fitems);
  /* optionally specify the _global_ numbers of the destination
     graph nodes of each forward item.
     this will define the reverse items (items2content) by
     insisting that these items be sorted in ascending order
     of their global number.
     reverse roots are not specified by this function, they
     are assumed to be the same as reverse items.
     one may only call this API or set_dest_idxs(), not both */
  void set_dest_globals(GOs fitems2ritem_globals);
  Dist invert() const;
  template <typename T>
  Read<T> exch(Read<T> data, Int width) const;
  template <typename T>
  Future<T> iexch(Read<T> data, Int width) const;
  template <typename T>
  Read<T> exch_reduce(Read<T> data, Int width, Omega_h_Op op) const;
  CommPtr parent_comm() const;
  CommPtr comm() const;
  LOs msgs2content() const;
  LOs content2msgs() const;
  LOs items2content() const;
  LOs items2msgs() const;
  LOs roots2items() const;
  Read<I32> msgs2ranks() const;
  Read<I32> items2ranks() const;
  LOs items2dest_idxs() const;
  Remotes items2dests() const;
  LO nitems() const;
  LO nroots() const;
  LO ndests() const;
  LO nsrcs() const;
  void change_comm(CommPtr new_comm);
  Remotes exch(Remotes data, Int width) const;

 private:
  void copy(Dist const& other);
  enum { F, R };
};

/*! \brief Creates a Dist object that can be re-used to synchronize variable-sized data per actor
   \param copies2owners Dist which sends data from each actor in this MPI rank to its owner
   \param copies2data A "fan" which, for each actor in this MPI rank, points to the start of its
                      variable-sized data in the array of all packed variable-sized data
   \details This function assumes and does not verify that the amount of variable-sized
            data per actor is consistent between MPI ranks.
 */
Dist create_dist_for_variable_sized(Dist copies2owners, LOs copies2data);

#define OMEGA_H_EXPL_INST_DECL(T)                                              \
  extern template Read<T> Dist::exch(Read<T> data, Int width) const;           \
  extern template Future<T> Dist::iexch(Read<T> data, Int width) const;        \
  extern template Read<T> Dist::exch_reduce<T>(                                \
      Read<T> data, Int width, Omega_h_Op op) const;
OMEGA_H_EXPL_INST_DECL(I8)
OMEGA_H_EXPL_INST_DECL(I32)
OMEGA_H_EXPL_INST_DECL(I64)
OMEGA_H_EXPL_INST_DECL(Real)
#undef OMEGA_H_EXPL_INST_DECL

}  // namespace Omega_h

#endif
