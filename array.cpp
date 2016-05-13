template <typename T>
Write<T>::Write(UInt size):
#ifdef USE_KOKKOS
  view_("omega_h", size)
#else
  ptr_(new T[size], std::default_delete<T[]>()),
  size_(size)
#endif
{
}

template <typename T>
static void fill(Write<T> a, T val) {
  auto f = LAMBDA(UInt i) {
    a[i] = val;
  };
  parallel_for(a.size(), f);
}

template <typename T>
Write<T>::Write(UInt size, T value):
  Write<T>(size)
{
  fill(*this, value);
}

template <typename T>
UInt Write<T>::size() const {
  return size_;
}

#ifdef USE_KOKKOS
template <typename T>
Kokkos::View<T*> Write<T>::view() const {
  return view_;
}
#endif

template <typename T>
Read<T>::Read() {
}

template <typename T>
Read<T>::Read(Write<T> write):
  write_(write) {
}

template <typename T>
Read<T>::Read(std::initializer_list<T> l):
  Read<T>(HostWrite<T>(l).write()) {
}

template <typename T>
UInt Read<T>::size() const {
  return write_.size();
}

#ifdef USE_KOKKOS
template <typename T>
Kokkos::View<const T*> Read<T>::view() const {
  return Kokkos::View<const T*>(write_.view());
}
#endif

template <typename T>
HostWrite<T>::HostWrite(UInt size):
  write_(size)
#ifdef USE_KOKKOS
  ,mirror_(Kokkos::create_mirror_view(write_.view()))
#endif
{
}

template <typename T>
HostWrite<T>::HostWrite(std::initializer_list<T> l):
  HostWrite<T>(l.size()) {
  UInt i;
  for (auto v : l)
    operator[](i++) = v;
}

template <typename T>
Write<T> HostWrite<T>::write() const {
#ifdef USE_KOKKOS
  Kokkos::deep_copy(write_.view(), mirror_);
#endif
  return write_;
}

template <typename T>
UInt HostWrite<T>::size() const {
  return write_.size();
}

template <typename T>
HostRead<T>::HostRead(Read<T> read):
  read_(read)
#ifdef USE_KOKKOS
  ,mirror_(Kokkos::create_mirror_view(read.view()))
#endif
{
#ifdef USE_KOKKOS
  Kokkos::deep_copy(mirror_, read_.view());
#endif
}

template <typename T>
UInt HostRead<T>::size() const {
  return read_.size();
}
