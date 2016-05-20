Mesh::Mesh():
  dim_(-1) {
  for (I8 i = 0; i <= 3; ++i)
    nents_[i] = -1;
}

void Mesh::set_dim(I8 dim) {
  CHECK(dim_ == -1);
  CHECK(dim >= 2);
  CHECK(dim <= 3);
  dim_ = dim;
}

void Mesh::set_verts(LO nverts) {
  nents_[VERT] = nverts;
}

void Mesh::set_ents(I8 dim, Adj down) {
  check_dim(dim);
  CHECK(!has_ents(dim));
  LOs hl2l = down.ab2b;
  CHECK(hl2l.size() % degrees[dim][dim - 1] == 0);
  nents_[dim] = hl2l.size() / degrees[dim][dim - 1];
  add_adj(dim, dim - 1, down);
  /* todo: replace with hilbert curve on set_coords() */
  add_tag_priv<GO>(dim, "global", 1, Read<GO>(nents(dim), 0, 1));
}

I8 Mesh::dim() const {
  return dim_;
}

LO Mesh::nents(I8 dim) const {
  check_dim2(dim);
  return nents_[dim];
}

template <typename T>
void Mesh::add_tag(I8 dim, std::string const& name, I8 ncomps, Read<T> data) {
  CHECK(!(dim == 0 && name == "coordinates"));
  add_tag_priv(dim, name, ncomps, data);
}

void Mesh::remove_tag(I8 dim, std::string const& name) {
  CHECK(!(dim == 0 && name == "coordinates"));
  remove_tag_priv(dim, name);
}

bool Mesh::has_tag(I8 dim, std::string const& name) const {
  check_dim2(dim);
  return tag_iter(dim, name) != tags_[dim].end();
}

template <typename T>
Tag<T> const& Mesh::get_tag(I8 dim, std::string const& name) const {
  check_dim2(dim);
  CHECK(has_tag(dim, name));
  return *(to<T>(tag_iter(dim, name)->get()));
}

I8 Mesh::count_tags(I8 dim) const {
  check_dim2(dim);
  return static_cast<I8>(tags_[dim].size());
}

TagBase const* Mesh::get_tag(I8 dim, I8 i) const {
  check_dim2(dim);
  CHECK(0 <= i);
  CHECK(i <= count_tags(dim));
  return tags_[dim][static_cast<std::size_t>(i)].get();
}

bool Mesh::has_ents(I8 dim) const {
  check_dim(dim);
  return nents_[dim] >= 0;
}

bool Mesh::has_adj(I8 from, I8 to) const {
  check_dim(from);
  check_dim(to);
  return bool(adjs_[from][to]);
}

Adj Mesh::get_adj(I8 from, I8 to) const {
  check_dim2(from);
  check_dim2(to);
  CHECK(has_adj(from, to));
  return *(adjs_[from][to]);
}

Adj Mesh::ask_adj(I8 from, I8 to) {
  check_dim2(from);
  check_dim2(to);
  if (has_adj(from, to))
    return get_adj(from,to);
  Adj derived = derive_adj(from, to);
  adjs_[from][to] = AdjPtr(new Adj(derived));
  return derived;
}

Read<GO> Mesh::ask_globals(I8 dim) {
  check_dim2(dim);
  if (!has_tag(dim, "global")) {
    /* todo: replace with hilbert curve */
    add_tag_priv<GO>(VERT, "global", 1, Read<GO>(nents(dim), 0, 1));
  }
  return get_tag<GO>(dim, "global").data();
}

struct HasName {
  std::string const& name_;
  HasName(std::string const& name):name_(name) {}
  bool operator()(std::shared_ptr<TagBase> const& a) {
    return a->name() == name_;
  }
};

template <typename T>
void Mesh::add_tag_priv(I8 dim, std::string const& name, I8 ncomps,
    Read<T> data) {
  check_dim2(dim);
  CHECK(data.size() == nents(dim) * ncomps);
  CHECK(!has_tag(dim, name));
  CHECK(tags_[dim].size() < static_cast<std::size_t>(INT8_MAX));
  tags_[dim].push_back(TagPtr(new Tag<T>(name, ncomps, data)));
}

void Mesh::remove_tag_priv(I8 dim, std::string const& name) {
  check_dim2(dim);
  CHECK(has_tag(dim, name));
  tags_[dim].erase(tag_iter(dim, name));
}

Mesh::TagIter Mesh::tag_iter(I8 dim, std::string const& name) {
  return std::find_if(begin(tags_[dim]), end(tags_[dim]), HasName(name));
}

Mesh::TagCIter Mesh::tag_iter(I8 dim, std::string const& name) const {
  return std::find_if(begin(tags_[dim]), end(tags_[dim]), HasName(name));
}

void Mesh::check_dim(I8 dim) const {
  CHECK(0 <= dim);
  CHECK(dim <= this->dim());
}

void Mesh::check_dim2(I8 dim) const {
  check_dim(dim);
  CHECK(has_ents(dim));
}

void Mesh::add_adj(I8 from, I8 to, Adj adj) {
  check_dim2(from);
  check_dim2(to);
  if (to < from) {
    CHECK(adj.a2ab.size() == 0);
    if (to == VERT)
      CHECK(adj.codes.size() == 0);
    CHECK(adj.ab2b.size() == nents(from) * degrees[from][to]);
  } else {
    if (from < to) {
      CHECK(adj.ab2b.size() == nents(to) * degrees[to][from]);
    }
    CHECK(adj.a2ab.size() == nents(from) + 1);
  }
  adjs_[from][to] = AdjPtr(new Adj(adj));
}

Adj Mesh::derive_adj(I8 from, I8 to) {
  check_dim(from);
  check_dim2(to);
  if (from < to) {
    Adj down = ask_adj(to, from);
    I8 nlows_per_high = degrees[to][from];
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
    /* todo: star code */
    NORETURN(Adj());
  }
}

#define INST_T(T) \
template Tag<T> const& Mesh::get_tag<T>( \
    I8 dim, std::string const& name) const; \
template void Mesh::add_tag(I8 dim, std::string const& name, \
    I8 ncomps, Read<T> data);
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T
