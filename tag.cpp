TagBase::TagBase(std::string const& name, I8 ncomps):
  name_(name),ncomps_(ncomps) {
}

TagBase::~TagBase() {
}

std::string const& TagBase::name() const {
  return name_;
}

I8 TagBase::ncomps() const {
  return ncomps_;
}

template <typename T>
bool is(TagBase const* t) {
  return nullptr != dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T> const* to(TagBase const* t) {
  CHECK(is<T>(t));
  return dynamic_cast<Tag<T> const*>(t);
}

template <typename T>
Tag<T>* to(TagBase* t) {
  CHECK(is<T>(t));
  return dynamic_cast<Tag<T>*>(t);
}

template <typename T>
Tag<T>::Tag(std::string const& name, I8 ncomps):
  TagBase(name, ncomps) {
}

template <typename T>
Read<T> Tag<T>::data() const {
  return data_;
}

template <typename T>
void Tag<T>::set_data(Read<T> data) {
  data_ = data;
}

#define INST_T(T) \
template bool is<T>(TagBase const* t); \
template Tag<T> const* to<T>(TagBase const* t); \
template class Tag<T>;
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T
