#include "Omega_h_dist.hpp"

#include "Omega_h_array_ops.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_map.hpp"
#include "Omega_h_int_scan.hpp"
#include "Omega_h_sort.hpp"

namespace Omega_h {

Dist::Dist() {}

Dist::Dist(Dist const& other) { copy(other); }

Dist& Dist::operator=(Dist const& other) {
  copy(other);
  return *this;
}

Dist::Dist(CommPtr comm_in, Remotes fitems2rroots, LO nrroots) {
  set_parent_comm(comm_in);
  set_dest_ranks(fitems2rroots.ranks);
  set_dest_idxs(fitems2rroots.idxs, nrroots);
}

void Dist::set_parent_comm(CommPtr parent_comm_in) {
  parent_comm_ = parent_comm_in;
}

void Dist::set_dest_ranks(Read<I32> items2ranks_in) {
  begin_code("Dist::set_dest_ranks");
  constexpr bool use_small_neighborhood_algorithm = true;
  if (use_small_neighborhood_algorithm) {
    Read<I32> msgs2ranks1;
    sort_small_range(
        items2ranks_in, &items2content_[F], &msgs2content_[F], &msgs2ranks1);
    comm_[F] = parent_comm_->graph(msgs2ranks1);
  } else {
    auto content2items = sort_by_keys(items2ranks_in);
    auto content2ranks = unmap(content2items, items2ranks_in, 1);
    Write<I8> jumps(content2ranks.size());
    auto mark_jumps = OMEGA_H_LAMBDA(LO i) {
      jumps[i] = (content2ranks[i] != content2ranks[i + 1]);
    };
    parallel_for(jumps.size() - 1, mark_jumps, "mark_jumps");
    if (jumps.size()) {
      jumps.set(jumps.size() - 1, 1);
    }
    auto tmp_content2msgs = offset_scan(Read<I8>(jumps));
    auto nmsgs = tmp_content2msgs.last();
    Write<I32> msgs2ranks_w(nmsgs);
    auto log_ranks = OMEGA_H_LAMBDA(LO i) {
      if (jumps[i]) {
        msgs2ranks_w[tmp_content2msgs[i]] = content2ranks[i];
      }
    };
    parallel_for(jumps.size(), log_ranks, "log_ranks");
    Write<LO> msgs2content_w(nmsgs + 1);
    msgs2content_w.set(0, 0);
    auto log_ends = OMEGA_H_LAMBDA(LO i) {
      if (jumps[i]) {
        msgs2content_w[tmp_content2msgs[i] + 1] = i + 1;
      }
    };
    parallel_for(jumps.size(), log_ends, "log_ends");
    items2content_[F] = invert_permutation(content2items);
    msgs2content_[F] = msgs2content_w;
    comm_[F] = parent_comm_->graph(msgs2ranks_w);
  }
  comm_[R] = comm_[F]->graph_inverse();
  auto fdegrees = get_degrees(msgs2content_[F]);
  auto rdegrees = comm_[F]->alltoall(fdegrees);
  msgs2content_[R] = offset_scan(rdegrees);
  end_code();
}

void Dist::set_dest_idxs(LOs fitems2rroots, LO nrroots) {
  begin_code("Dist::set_dest_idxs");
  auto rcontent2rroots = exch(fitems2rroots, 1);
  auto rroots2rcontent = invert_map_by_atomics(rcontent2rroots, nrroots);
  roots2items_[R] = rroots2rcontent.a2ab;
  items2content_[R] = rroots2rcontent.ab2b;
  end_code();
}

void Dist::set_dest_globals(GOs fitems2ritem_globals) {
  begin_code("Dist::set_dest_globals");
  auto rcontent2ritem_globals = exch(fitems2ritem_globals, 1);
  items2content_[R] = sort_by_keys(rcontent2ritem_globals);
  roots2items_[R] = LOs();
  end_code();
}

void Dist::set_roots2items(LOs froots2fitems) {
  roots2items_[F] = froots2fitems;
}

Dist Dist::invert() const {
  Dist out;
  out.parent_comm_ = parent_comm_;
  for (Int i = 0; i < 2; ++i) {
    out.roots2items_[i] = roots2items_[1 - i];
    out.items2content_[i] = items2content_[1 - i];
    out.msgs2content_[i] = msgs2content_[1 - i];
    out.comm_[i] = comm_[1 - i];
  }
  return out;
}

template <typename T>
Read<T> Dist::exch(Read<T> data, Int width) const {
  ScopedTimer exch_timer("Dist::exch");
  if (roots2items_[F].exists()) {
    data = expand(data, roots2items_[F], width);
  }
  if (items2content_[F].exists()) {
    data = permute(data, items2content_[F], width);
  }
  data = comm_[F]->alltoallv(data, msgs2content_[F], msgs2content_[R], width);
  if (items2content_[R].exists()) {
    data = unmap(items2content_[R], data, width);
  }
  return data;
}

template <typename T>
Read<T> Dist::exch_reduce(Read<T> data, Int width, Omega_h_Op op) const {
  Read<T> item_data = exch(data, width);
  return fan_reduce(roots2items_[R], item_data, width, op);
}

CommPtr Dist::parent_comm() const { return parent_comm_; }

CommPtr Dist::comm() const { return comm_[F]; }

LOs Dist::msgs2content() const { return msgs2content_[F]; }

LOs Dist::content2msgs() const { return invert_fan(msgs2content_[F]); }

LOs Dist::items2msgs() const {
  return unmap(items2content_[F], content2msgs(), 1);
}

LOs Dist::roots2items() const { return roots2items_[F]; }

Read<I32> Dist::msgs2ranks() const { return comm_[F]->destinations(); }

Read<I32> Dist::items2ranks() const {
  return compound_maps(items2msgs(), msgs2ranks());
}

LOs Dist::items2dest_idxs() const {
  auto inverse = invert();
  return inverse.exch(LOs(ndests(), 0, 1), 1);
}

Remotes Dist::items2dests() const {
  return Remotes(items2ranks(), items2dest_idxs());
}

LO Dist::nitems() const { return msgs2content_[F].last(); }

LO Dist::nroots() const { return roots2items_[F].size() - 1; }

LO Dist::nsrcs() const {
  if (roots2items_[F].exists()) return nroots();
  return nitems();
}

LO Dist::ndests() const { return invert().nsrcs(); }

/* this is the key algorithm for moving from one communicator
   to another. essentially, we have to map from old ranks to
   new ranks, and rebuild graph communicators as well */
void Dist::change_comm(CommPtr new_comm) {
  // gather the new ranks of our neighbors
  auto new_sources = comm_[F]->allgather(new_comm->rank());
  auto new_destinations = comm_[R]->allgather(new_comm->rank());
  // rebuild graph communicators from these new neighbor lists
  comm_[F] = new_comm->graph_adjacent(new_sources, new_destinations);
  comm_[R] = comm_[F]->graph_inverse();
  // replace parent_comm_
  parent_comm_ = new_comm;
  // thats it! since all rank information is queried from graph comms
}

Remotes Dist::exch(Remotes data, Int width) const {
  auto ranks = exch(data.ranks, width);
  auto idxs = exch(data.idxs, width);
  return Remotes(ranks, idxs);
}

void Dist::copy(Dist const& other) {
  parent_comm_ = other.parent_comm_;
  for (Int i = 0; i < 2; ++i) {
    roots2items_[i] = other.roots2items_[i];
    items2content_[i] = other.items2content_[i];
    msgs2content_[i] = other.msgs2content_[i];
    comm_[i] = other.comm_[i];
  }
}

#define INST_T(T)                                                              \
  template Read<T> Dist::exch(Read<T> data, Int width) const;                  \
  template Read<T> Dist::exch_reduce(Read<T> data, Int width, Omega_h_Op op)   \
      const;
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T

}  // end namespace Omega_h
