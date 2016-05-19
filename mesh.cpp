Mesh::Mesh():
  dim_(-1) {
  for (I8 i = 0; i <= 3; ++i)
    nents_[i] = -1;
}

Mesh::Mesh(I8 dim, LO nverts):
  dim_(dim) {
  CHECK(dim > 0);
  CHECK(dim <= 3);
  nents_[VERT] = nverts;
  for (I8 i = 1; i <= dim; ++i)
    nents_[i] = -1;
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
  check_dim2(dim);
  CHECK(!has_tag(dim, name));
  CHECK(tags_[dim].size() < static_cast<std::size_t>(INT8_MAX));
  tags_[dim].push_back(TagPtr(new Tag<T>(name, ncomps, data)));
}

void Mesh::remove_tag(I8 dim, std::string const& name) {
  check_dim2(dim);
  CHECK(has_tag(dim, name));
  tags_[dim].erase(tag_iter(dim, name));
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

bool Mesh::has_dim(I8 dim) const {
  check_dim(dim);
  return nents_[dim] != -1;
}

void Mesh::add_adj(I8 from, I8 to, Adj adj) {
  check_dim(from);
  check_dim2(to);
  CHECK(!has_dim(from));
  CHECK(from > 0);
  CHECK((to == VERT) || (to == (from - 1)));
  adjs_[from][to] = AdjPtr(new Adj(adj));
}

bool Mesh::has_adj(I8 from, I8 to) const {
  check_dim(from);
  check_dim(to);
  return bool(adjs_[DIMS][DIMS]);
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

struct HasName {
  std::string const& name_;
  HasName(std::string const& name):name_(name) {}
  bool operator()(std::shared_ptr<TagBase> const& a) {
    return a->name() == name_;
  }
};

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
  CHECK(has_dim(dim));
}

Adj Mesh::derive_adj(I8 from, I8 to) {
  check_dim(from);
  check_dim2(to);
  /* todo... */
  return Adj();
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
