Dist::Dist() {
}

Dist::Dist(CommPtr comm, Remotes fitems2rroots, LO nrroots) {
  set_parent_comm(comm);
  set_dest_ranks(fitems2rroots.ranks);
  set_dest_idxs(fitems2rroots.idxs, nrroots);
}

void Dist::set_parent_comm(CommPtr parent_comm) {
  parent_comm_ = parent_comm;
}

void Dist::set_dest_ranks(Read<I32> items2ranks) {
  auto content2items = sort_by_keys(items2ranks);
  auto content2ranks = unmap(content2items, items2ranks, 1);
  Write<I8> jumps(content2ranks.size());
  auto mark_jumps = LAMBDA(LO i) {
    jumps[i] = (content2ranks[i] != content2ranks[i + 1]);
  };
  parallel_for(jumps.size() - 1, mark_jumps);
  if (jumps.size()) {
    jumps.set(jumps.size() - 1, 1);
  }
  auto content2msgs = offset_scan(Read<I8>(jumps));
  auto nmsgs = content2msgs.last();
  Write<I32> msgs2ranks(nmsgs);
  auto log_ranks = LAMBDA(LO i) {
    if (jumps[i]) {
      msgs2ranks[content2msgs[i]] = content2ranks[i];
    }
  };
  parallel_for(jumps.size(), log_ranks);
  Write<LO> msgs2content(nmsgs + 1);
  msgs2content.set(0, 0);
  auto log_ends = LAMBDA(LO i) {
    if (jumps[i]) {
      msgs2content[content2msgs[i] + 1] = i + 1;
    }
  };
  parallel_for(jumps.size(), log_ends);
  items2content_[F] = invert_permutation(content2items);
  msgs2content_[F] = msgs2content;
  msgs2ranks_[F] = msgs2ranks;
  comm_[F] = parent_comm_->graph(msgs2ranks);
  comm_[R] = parent_comm_->graph_adjacent(
      comm_[F]->destinations(), comm_[F]->sources());
  msgs2ranks_[R] = comm_[F]->sources();
  auto fdegrees = get_degrees(msgs2content_[F]);
  auto rdegrees = comm_[F]->alltoall(fdegrees);
  msgs2content_[R] = offset_scan(rdegrees);
}

void Dist::set_dest_idxs(LOs fitems2rroots, LO nrroots) {
  auto rcontent2rroots = exch(fitems2rroots, 1);
  map::invert_by_sorting(rcontent2rroots, nrroots,
      roots2items_[R], items2content_[R]);
}

void Dist::set_roots2items(LOs froots2fitems) {
  CHECK(froots2fitems.last() == nitems());
  roots2items_[F] = froots2fitems;
}

Dist Dist::invert() const {
  Dist out;
  out.parent_comm_ = parent_comm_;
  for (Int i = 0; i < 2; ++i) {
    out.roots2items_[i] = roots2items_[1 - i];
    out.items2content_[i] = items2content_[1 - i];
    out.msgs2content_[i] = msgs2content_[1 - i];
    out.msgs2ranks_[i] = msgs2ranks_[1 - i];
    out.comm_[i] = comm_[1 - i];
  }
  return out;
}

template <typename T>
Read<T> Dist::exch(Read<T> data, Int width) const {
  if (roots2items_[F].exists()) {
    data = expand(data, roots2items_[F], width);
  }
  if (items2content_[F].exists()) {
    data = permute(data, items2content_[F], width);
  }
  auto sendcounts = multiply_each_by(width,
      get_degrees(msgs2content_[F]));
  auto recvcounts = multiply_each_by(width,
      get_degrees(msgs2content_[R]));
  auto sdispls = offset_scan(sendcounts);
  auto rdispls = offset_scan(recvcounts);
  data = comm_[F]->alltoallv(data,
      sendcounts, sdispls,
      recvcounts, rdispls);
  if (items2content_[R].exists()) {
    data = unmap(items2content_[R], data, width);
  }
  return data;
}

template Read<I8> Dist::exch(Read<I8> data, Int width) const;
template Read<I32> Dist::exch(Read<I32> data, Int width) const;
template Read<I64> Dist::exch(Read<I64> data, Int width) const;
template Read<Real> Dist::exch(Read<Real> data, Int width) const;

CommPtr Dist::parent_comm() const {
  return parent_comm_;
}

CommPtr Dist::comm() const {
  return comm_[F];
}

LOs Dist::content2msgs() const {
  return invert_fan(msgs2content_[F]);
}

LOs Dist::items2msgs() const {
  return unmap(items2content_[F], content2msgs(), 1);
}

LOs Dist::roots2items() const {
  return roots2items_[F];
}

Read<I32> Dist::msgs2ranks() const {
  return msgs2ranks_[F];
}

Read<I32> Dist::items2ranks() const {
  return compound_maps(items2msgs(), msgs2ranks());
}

LO Dist::nitems() const {
  return items2content_[F].size();
}

LO Dist::nroots() const {
  return roots2items_[F].size() - 1;
}

void Dist::change_comm(CommPtr new_comm) {
  auto new_sources = comm_[F]->allgather(new_comm->rank());
  auto new_destinations = comm_[R]->allgather(new_comm->rank());
  comm_[F] = new_comm->graph_adjacent(new_sources, new_destinations);
  comm_[R] = new_comm->graph_adjacent(new_destinations, new_sources);
  parent_comm_ = new_comm;
}
