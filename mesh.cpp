Mesh::Mesh():
  dim_(-1),
  partition_(-1) {
  for (Int i = 0; i <= 3; ++i)
    nents_[i] = -1;
  partition_ = ELEMENT_BASED;
}

void Mesh::set_comm(CommPtr new_comm) {
  auto rank_had_comm = bool(comm_);
  auto nnew_had_comm = new_comm->allreduce(I32(rank_had_comm), SUM);
  if (0 < nnew_had_comm && nnew_had_comm < new_comm->size()) {
    //partitioning out from small sub-communicator to larger one
    if (!rank_had_comm) {
      //temporarily set the uninit ranks to Comm::self()
      comm_ = Comm::self();
    }
    bcast_mesh(*this, new_comm, rank_had_comm);
  }
  if (0 < nnew_had_comm) {
    for (Int d = 0; d <= dim(); ++d) {
      ask_globals(d);
      auto dist = ask_dist(d);
      dist.change_comm(new_comm);
      owners_[d].ranks = dist.items2ranks();
    }
  }
  comm_ = new_comm;
}

void Mesh::set_dim(Int dim) {
  CHECK(dim_ == -1);
  CHECK(dim >= 2);
  CHECK(dim <= 3);
  dim_ = dim;
}

void Mesh::set_verts(LO nverts) {
  nents_[VERT] = nverts;
}

void Mesh::set_ents(Int dim, Adj down) {
  check_dim(dim);
  CHECK(!has_ents(dim));
  LOs hl2l = down.ab2b;
  CHECK(hl2l.size() % simplex_degrees[dim][dim - 1] == 0);
  nents_[dim] = hl2l.size() / simplex_degrees[dim][dim - 1];
  add_adj(dim, dim - 1, down);
}

CommPtr Mesh::comm() const {
  return comm_;
}

Int Mesh::dim() const {
  return dim_;
}

LO Mesh::nents(Int dim) const {
  check_dim2(dim);
  return nents_[dim];
}

LO Mesh::nelems() const{
  return nents(dim());
}

LO Mesh::nverts() const {
  return nents(VERT);
}

LO Mesh::nedges() const {
  return nents(EDGE);
}

template <typename T>
void Mesh::add_tag(Int dim, std::string const& name, Int ncomps,
    Xfer xfer) {
  check_dim2(dim);
  CHECK(!has_tag(dim, name));
  CHECK(ncomps >= 0);
  CHECK(tags_[dim].size() < static_cast<std::size_t>(INT8_MAX));
  tags_[dim].push_back(TagPtr(new Tag<T>(name, ncomps, xfer)));
}

template <typename T>
void Mesh::add_tag(Int dim, std::string const& name, Int ncomps,
    Xfer xfer, Read<T> array) {
  add_tag<T>(dim, name, ncomps, xfer);
  set_tag<T>(dim, name, array);
}

template <typename T>
void Mesh::set_tag(Int dim, std::string const& name, Read<T> array) {
  CHECK(has_tag(dim, name));
  Tag<T>* tag = to<T>(tag_iter(dim, name)->get());
  CHECK(array.size() == nents(dim) * tag->ncomps());
  tag->set_array(array);
  react_to_set_tag(dim, name);
}

void Mesh::react_to_set_tag(Int dim, std::string const& name) {
  /* hardcoded cache invalidations */
  if ((dim == VERT) && ((name == "coordinates") ||
                        (name == "size") ||
                        (name == "metric"))) {
    if (has_tag(EDGE, "length")) remove_tag(EDGE, "length");
    if (has_tag(this->dim(), "quality")) remove_tag(this->dim(), "quality");
  }
}

template <typename T>
Tag<T> const& Mesh::get_tag(Int dim, std::string const& name) const {
  check_dim2(dim);
  if (!has_tag(dim, name))
    fail("expected tag %s on dimension %d\n", name.c_str(), dim);
  return *(to<T>(tag_iter(dim, name)->get()));
}

template <typename T>
Read<T> Mesh::get_array(Int dim, std::string const& name) const {
  return get_tag<T>(dim, name).array();
}

void Mesh::remove_tag(Int dim, std::string const& name) {
  check_dim2(dim);
  CHECK(has_tag(dim, name));
  tags_[dim].erase(tag_iter(dim, name));
}

bool Mesh::has_tag(Int dim, std::string const& name) const {
  check_dim(dim);
  if (!has_ents(dim))
    return false;
  return tag_iter(dim, name) != tags_[dim].end();
}

Int Mesh::ntags(Int dim) const {
  check_dim2(dim);
  return static_cast<Int>(tags_[dim].size());
}

TagBase const* Mesh::get_tag(Int dim, Int i) const {
  check_dim2(dim);
  CHECK(0 <= i);
  CHECK(i <= ntags(dim));
  return tags_[dim][static_cast<std::size_t>(i)].get();
}

bool Mesh::has_ents(Int dim) const {
  check_dim(dim);
  return nents_[dim] >= 0;
}

bool Mesh::has_adj(Int from, Int to) const {
  check_dim(from);
  check_dim(to);
  return bool(adjs_[from][to]);
}

Adj Mesh::get_adj(Int from, Int to) const {
  check_dim2(from);
  check_dim2(to);
  CHECK(has_adj(from, to));
  return *(adjs_[from][to]);
}

Adj Mesh::ask_down(Int from, Int to) {
  CHECK(to < from);
  return ask_adj(from, to);
}

LOs Mesh::ask_verts_of(Int dim) {
  return ask_adj(dim, VERT).ab2b;
}

Adj Mesh::ask_up(Int from, Int to) {
  CHECK(from < to);
  return ask_adj(from, to);
}

Graph Mesh::ask_star(Int dim) {
  CHECK(dim < this->dim());
  return ask_adj(dim, dim);
}

Graph Mesh::ask_dual() {
  return ask_adj(dim(), dim());
}

struct HasName {
  std::string const& name_;
  HasName(std::string const& name):name_(name) {}
  bool operator()(std::shared_ptr<TagBase> const& a) {
    return a->name() == name_;
  }
};

Mesh::TagIter Mesh::tag_iter(Int dim, std::string const& name) {
  return std::find_if(tags_[dim].begin(), tags_[dim].end(), HasName(name));
}

Mesh::TagCIter Mesh::tag_iter(Int dim, std::string const& name) const {
  return std::find_if(tags_[dim].begin(), tags_[dim].end(), HasName(name));
}

void Mesh::check_dim(Int dim) const {
  CHECK(0 <= dim);
  CHECK(dim <= this->dim());
}

void Mesh::check_dim2(Int dim) const {
  check_dim(dim);
  CHECK(has_ents(dim));
}

void Mesh::add_adj(Int from, Int to, Adj adj) {
  check_dim2(from);
  check_dim(to);
  if (to < from) {
    CHECK(!adj.a2ab.exists());
    if (to == VERT)
      CHECK(!adj.codes.exists());
    CHECK(adj.ab2b.size() == nents(from) * simplex_degrees[from][to]);
  } else {
    if (from < to) {
      CHECK(adj.ab2b.size() == nents(to) * simplex_degrees[to][from]);
    }
    CHECK(adj.a2ab.size() == nents(from) + 1);
  }
  adjs_[from][to] = AdjPtr(new Adj(adj));
}

Adj Mesh::derive_adj(Int from, Int to) {
  check_dim(from);
  check_dim2(to);
  if (from < to) {
    Adj down = ask_adj(to, from);
    Int nlows_per_high = simplex_degrees[to][from];
    LO nlows = nents(from);
    Read<GO> high_globals = ask_globals(to);
    Adj up = invert(down, nlows_per_high, nlows, high_globals);
    return up;
  } else if (to < from) {
    CHECK(to + 1 < from);
    Adj h2m = ask_adj(from, to + 1);
    Adj m2l = ask_adj(to + 1, to);
    Adj h2l = transit(h2m, m2l, from, to);
    return h2l;
  } else {
    if (from == dim() && to == dim()) {
      return elements_across_sides(dim(),
          ask_adj(dim(), dim() - 1), ask_adj(dim() - 1, dim()),
          mark_exposed_sides(*this));
    }
    if (from == VERT && to == VERT) {
      return verts_across_edges(ask_adj(EDGE,VERT), ask_adj(VERT,EDGE));
    }
    if (from == EDGE && to == EDGE) {
      CHECK(dim() >= 2);
      Graph g = edges_across_tris(ask_adj(TRI,EDGE),ask_adj(EDGE,TRI));
      if (dim() == 3) {
        g = add_edges(g, edges_across_tets(
              ask_adj(TET,EDGE), ask_adj(EDGE,TET)));
      }
      return g;
    }
  }
  fail("can't derive adjacency from %s to %s\n",
      plural_names[from], plural_names[to]);
}

Adj Mesh::ask_adj(Int from, Int to) {
  check_dim2(from);
  check_dim2(to);
  if (has_adj(from, to)) {
    return get_adj(from,to);
  }
  Adj derived = derive_adj(from, to);
  adjs_[from][to] = AdjPtr(new Adj(derived));
  return derived;
}

void Mesh::add_coords(Reals array) {
  add_tag<Real>(0, "coordinates", dim(), OSH_LINEAR_INTERP, array);
}

Reals Mesh::coords() const {
  return get_array<Real>(0, "coordinates");
}

void Mesh::set_coords(Reals array) {
  set_tag<Real>(0, "coordinates", array);
}

Read<GO> Mesh::ask_globals(Int dim) {
  if (!has_tag(dim, "global")) {
    CHECK(comm_->size() == 1);
    add_tag(dim, "global", 1, OSH_GLOBAL, Read<GO>(nents(dim), 0, 1));
  }
  return get_array<GO>(dim, "global");
}

void Mesh::forget_globals() {
  for (Int d = 0; d <= dim(); ++d)
    if (has_tag(d, "global"))
      remove_tag(d, "global");
}

Reals Mesh::ask_edge_lengths() {
  if (!has_tag(EDGE, "length")) {
    auto lengths = measure_edges(*this);
    add_tag(EDGE, "length", 1, OSH_LENGTH, lengths);
  }
  return get_array<Real>(EDGE, "length");
}

Reals Mesh::ask_qualities() {
  if (!has_tag(dim(), "quality")) {
    auto qualities = measure_qualities(*this);
    add_tag(dim(), "quality", 1, OSH_QUALITY, qualities);
  }
  return get_array<Real>(dim(), "quality");
}

void Mesh::set_owners(Int dim, Remotes owners) {
  owners_[dim] = owners;
  dists_[dim] = DistPtr();
}

Remotes Mesh::ask_owners(Int dim) {
  if (!owners_[dim].ranks.exists() ||
      !owners_[dim].idxs.exists()) {
    CHECK(comm_->size() == 1);
    owners_[dim] = Remotes(Read<I32>(nents(dim), comm_->rank()),
        LOs(nents(dim), 0, 1));
  }
  return owners_[dim];
}

Read<I8> Mesh::owned(Int dim) {
  auto e2rank = ask_owners(dim).ranks;
  return mark_equal(e2rank, comm()->rank());
}

Dist Mesh::ask_dist(Int dim) {
  if (!dists_[dim]) {
    auto owners = ask_owners(dim);
    CHECK(owners.ranks.exists());
    CHECK(owners.idxs.exists());
    dists_[dim] = DistPtr(new Dist(comm_, owners, nents(dim)));
  }
  return *(dists_[dim]);
}

Partition Mesh::partition() const {
  CHECK(partition_ != -1);
  return static_cast<Partition>(partition_);
}

void Mesh::set_partition(Partition partition) {
  if (partition_ == -1) {
    partition_ = partition;
    return;
  }
  if (partition_ == partition) {
    return;
  }
  if (partition_ != ELEMENT_BASED) {
    partition_by_elems(*this);
    partition_ = ELEMENT_BASED;
  }
  if (partition == GHOSTED) {
    ghost_mesh(*this);
  } else if (partition == VERTEX_BASED) {
    partition_by_verts(*this);
  }
  partition_ = partition;
}

void Mesh::migrate(Remotes new_elems2old_owners) {
  migrate_mesh(*this, new_elems2old_owners);
}

void Mesh::reorder() {
  reorder_by_hilbert(*this);
}

void Mesh::balance() {
  set_partition(ELEMENT_BASED);
  inertia::Rib hints;
  if (rib_hints_)
    hints = *rib_hints_;
  auto ecoords = average_field(*this, dim(), LOs(nelems(), 0, 1), dim(),
      coords());
  if (dim() == 2)
    ecoords = vectors_2d_to_3d(ecoords);
  auto masses = Reals(nelems(), 1);
  auto owners = ask_owners(dim());
  auto total = comm_->allreduce(GO(nelems()), SUM);
  auto avg = Real(total) / Real(comm_->size());
  hints = recursively_bisect(comm(), ecoords, masses, owners,
      2.0 / avg, hints);
  rib_hints_ = RibPtr(new inertia::Rib(hints));
  migrate(owners);
}

Graph Mesh::ask_graph(Int from, Int to) {
  if (to > from) {
    return ask_up(from, to);
  }
  if (to < from) {
    auto down = ask_down(from, to);
    auto a2ab = LOs(nents(from) + 1, 0, simplex_degrees[from][to]);
    return Graph(a2ab, down.ab2b);
  }
  CHECK(from == to);
  return Graph(LOs(nents(to) + 1, 0, 1), LOs(nents(to), 0, 1));
}

template <typename T>
Read<T> Mesh::sync_array(Int ent_dim, Read<T> a, Int width) {
  return ask_dist(ent_dim).invert().exch(a, width);
}

#define INST_T(T) \
template Tag<T> const& Mesh::get_tag<T>( \
    Int dim, std::string const& name) const; \
template Read<T> Mesh::get_array<T>( \
    Int dim, std::string const& name) const; \
template void Mesh::add_tag<T>(Int dim, std::string const& name, Int ncomps, \
    Xfer xfer); \
template void Mesh::add_tag<T>(Int dim, std::string const& name, Int ncomps, \
    Xfer xfer, Read<T> array); \
template void Mesh::set_tag(Int dim, std::string const& name, Read<T> array); \
template Read<T> Mesh::sync_array(Int ent_dim, Read<T> a, Int width);
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T
