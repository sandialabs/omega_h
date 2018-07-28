#include "Omega_h_comm.hpp"

#include <string>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_scan.hpp"

#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
#include "Omega_h_library.hpp"
#include "Omega_h_for.hpp"
#endif

#include "Omega_h_map.hpp"
#include "Omega_h_stack.hpp"

namespace Omega_h {

#ifdef OMEGA_H_USE_MPI
#define CALL(f) OMEGA_H_CHECK(MPI_SUCCESS == (f))
#endif

Comm::Comm() {
#ifdef OMEGA_H_USE_MPI
  impl_ = MPI_COMM_NULL;
#endif
  library_ = nullptr;
}

#ifdef OMEGA_H_USE_MPI
Comm::Comm(Library* library_in, MPI_Comm impl_in)
    : impl_(impl_in), library_(library_in) {
  int topo_type;
  CALL(MPI_Topo_test(impl_in, &topo_type));
  if (topo_type == MPI_DIST_GRAPH) {
    int nin, nout, is_weighted;
    CALL(MPI_Dist_graph_neighbors_count(impl_in, &nin, &nout, &is_weighted));
    HostWrite<I32> sources_w(nin);
    HostWrite<I32> destinations_w(nout);
    CALL(MPI_Dist_graph_neighbors(impl_in, nin, nonnull(sources_w.data()),
        OMEGA_H_MPI_UNWEIGHTED, nout, nonnull(destinations_w.data()),
        OMEGA_H_MPI_UNWEIGHTED));
    srcs_ = sources_w.write();
    dsts_ = destinations_w.write();
    self_src_ = find_last(srcs_, rank());
    self_dst_ = find_last(dsts_, rank());
    host_srcs_ = HostRead<I32>(srcs_);
    host_dsts_ = HostRead<I32>(dsts_);
  }
}
#else
Comm::Comm(Library* library_in, bool is_graph, bool sends_to_self)
    : library_(library_in) {
  if (is_graph) {
    if (sends_to_self) {
      srcs_ = Read<LO>({0});
      self_src_ = self_dst_ = 0;
    } else {
      srcs_ = Read<LO>({});
      self_src_ = self_dst_ = -1;
    }
    dsts_ = srcs_;
    host_srcs_ = HostRead<I32>(srcs_);
    host_dsts_ = HostRead<I32>(dsts_);
  } else {
    OMEGA_H_CHECK(!sends_to_self);
  }
}
#endif

Comm::~Comm() {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Comm_free(&impl_));
#endif
}

Library* Comm::library() const { return library_; }

I32 Comm::rank() const {
#ifdef OMEGA_H_USE_MPI
  I32 r;
  CALL(MPI_Comm_rank(impl_, &r));
  return r;
#else
  return 0;
#endif
}

I32 Comm::size() const {
#ifdef OMEGA_H_USE_MPI
  I32 s;
  CALL(MPI_Comm_size(impl_, &s));
  return s;
#else
  return 1;
#endif
}

CommPtr Comm::dup() const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  CALL(MPI_Comm_dup(impl_, &impl2));
  return CommPtr(new Comm(library_, impl2));
#else
  return CommPtr(
      new Comm(library_, srcs_.exists(), srcs_.exists() && srcs_.size() == 1));
#endif
}

CommPtr Comm::split(I32 color, I32 key) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  CALL(MPI_Comm_split(impl_, color, key, &impl2));
  return CommPtr(new Comm(library_, impl2));
#else
  (void)color;
  (void)key;
  return CommPtr(new Comm(library_, false, false));
#endif
}

CommPtr Comm::graph(Read<I32> dsts) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  int n = 1;
  int source[1] = {rank()};
  int degree[1] = {dsts.size()};
  HostRead<I32> h_destinations(dsts);
  int reorder = 0;
  CALL(MPI_Dist_graph_create(impl_, n, source, degree,
      nonnull(h_destinations.data()), OMEGA_H_MPI_UNWEIGHTED, MPI_INFO_NULL,
      reorder, &impl2));
  return CommPtr(new Comm(library_, impl2));
#else
  return CommPtr(new Comm(library_, true, dsts.size() == 1));
#endif
}

CommPtr Comm::graph_adjacent(Read<I32> srcs, Read<I32> dsts) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  HostRead<I32> h_sources(srcs);
  HostRead<I32> h_destinations(dsts);
  int reorder = 0;
  CALL(MPI_Dist_graph_create_adjacent(impl_, h_sources.size(),
      nonnull(h_sources.data()), OMEGA_H_MPI_UNWEIGHTED, h_destinations.size(),
      nonnull(h_destinations.data()), OMEGA_H_MPI_UNWEIGHTED, MPI_INFO_NULL,
      reorder, &impl2));
  return CommPtr(new Comm(library_, impl2));
#else
  OMEGA_H_CHECK(srcs == dsts);
  return CommPtr(new Comm(library_, true, dsts.size() == 1));
#endif
}

CommPtr Comm::graph_inverse() const {
  return graph_adjacent(destinations(), sources());
}

Read<I32> Comm::sources() const { return srcs_; }

Read<I32> Comm::destinations() const { return dsts_; }

template <typename T>
T Comm::allreduce(T x, Omega_h_Op op) const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Allreduce(
      MPI_IN_PLACE, &x, 1, MpiTraits<T>::datatype(), mpi_op(op), impl_));
#else
  (void)op;
#endif
  return x;
}

bool Comm::reduce_or(bool x) const {
  I8 y = x;
  y = allreduce(y, OMEGA_H_MAX);
  return static_cast<bool>(y);
}

bool Comm::reduce_and(bool x) const {
  I8 y = x;
  y = allreduce(y, OMEGA_H_MIN);
  return static_cast<bool>(y);
}

#ifdef OMEGA_H_USE_MPI
static void mpi_add_int128(void* a, void* b, int*, MPI_Datatype*) {
  Int128* a2 = static_cast<Int128*>(a);
  Int128* b2 = static_cast<Int128*>(b);
  *b2 = *b2 + *a2;
}
#endif

Int128 Comm::add_int128(Int128 x) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Op op;
  int commute = true;
  CALL(MPI_Op_create(mpi_add_int128, commute, &op));
  CALL(MPI_Allreduce(MPI_IN_PLACE, &x, sizeof(Int128), MPI_PACKED, op, impl_));
  CALL(MPI_Op_free(&op));
#endif
  return x;
}

template <typename T>
T Comm::exscan(T x, Omega_h_Op op) const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Exscan(
      MPI_IN_PLACE, &x, 1, MpiTraits<T>::datatype(), mpi_op(op), impl_));
  if (rank() == 0) x = 0;
  return x;
#else
  (void)op;
  (void)x;
  return 0;
#endif
}

template <typename T>
void Comm::bcast(T& x) const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Bcast(&x, 1, MpiTraits<T>::datatype(), 0, impl_));
#else
  (void)x;
#endif
}

void Comm::bcast_string(std::string& s) const {
#ifdef OMEGA_H_USE_MPI
  I32 len = static_cast<I32>(s.length());
  bcast(len);
  s.resize(static_cast<std::size_t>(len));
  CALL(MPI_Bcast(&s[0], len, MPI_CHAR, 0, impl_));
#else
  (void)s;
#endif
}

#ifdef OMEGA_H_USE_MPI

/* custom implementation of MPI_Neighbor_allgather
 * in the case that we are using an MPI older
 * than version 3.0
 */

static int Neighbor_allgather(HostRead<I32> sources, HostRead<I32> destinations,
    const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
    int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
#if MPI_VERSION < 3
  static int const tag = 42;
  int indegree, outdegree;
  indegree = sources.size();
  outdegree = destinations.size();
  int recvwidth;
  CALL(MPI_Type_size(sendtype, &recvwidth));
  MPI_Request* recvreqs = new MPI_Request[indegree];
  MPI_Request* sendreqs = new MPI_Request[outdegree];
  for (int i = 0; i < indegree; ++i)
    CALL(MPI_Irecv(static_cast<char*>(recvbuf) + i * recvwidth, recvcount,
        recvtype, sources[i], tag, comm, recvreqs + i));
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(sendbuf, sendcount, sendtype, destinations[i], tag, comm,
        sendreqs + i));
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  delete[] sendreqs;
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  delete[] recvreqs;
  return MPI_SUCCESS;
#else
  (void)sources;
  (void)destinations;
  return MPI_Neighbor_allgather(
      sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
#endif  // end if MPI_VERSION < 3
}

/* custom implementation of MPI_Neighbor_alltoall
 * in the case that we are using an MPI older
 * than version 3.0
 */

static int Neighbor_alltoall(HostRead<I32> sources, HostRead<I32> destinations,
    const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
    int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
#if MPI_VERSION < 3
  static int const tag = 42;
  int indegree, outdegree;
  indegree = sources.size();
  outdegree = destinations.size();
  int sendwidth;
  CALL(MPI_Type_size(sendtype, &sendwidth));
  int recvwidth;
  CALL(MPI_Type_size(sendtype, &recvwidth));
  MPI_Request* recvreqs = new MPI_Request[indegree];
  MPI_Request* sendreqs = new MPI_Request[outdegree];
  for (int i = 0; i < indegree; ++i)
    CALL(MPI_Irecv(static_cast<char*>(recvbuf) + i * recvwidth, recvcount,
        recvtype, sources[i], tag, comm, recvreqs + i));
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(static_cast<char const*>(sendbuf) + i * sendwidth, sendcount,
        sendtype, destinations[i], tag, comm, sendreqs + i));
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  delete[] sendreqs;
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  delete[] recvreqs;
  return MPI_SUCCESS;
#else
  (void)sources;
  (void)destinations;
  return MPI_Neighbor_alltoall(
      sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
#endif  // end if MPI_VERSION < 3
}

/* custom implementation of MPI_Neighbor_alltoallv
 * this used to be here as a workaround, but now allows us to precompute
 * fewer things
 */

static int Neighbor_alltoallv(HostRead<I32> sources, HostRead<I32> destinations,
    int width, const void* sendbuf, const int sdispls[], MPI_Datatype sendtype,
    void* recvbuf, const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm) {
  begin_code("Neighbor_alltoallv");
  int const tag = 42;
  int indegree, outdegree;
  indegree = sources.size();
  outdegree = destinations.size();
  int sendwidth;
  CALL(MPI_Type_size(sendtype, &sendwidth));
  int recvwidth;
  CALL(MPI_Type_size(sendtype, &recvwidth));
  MPI_Request* recvreqs = new MPI_Request[indegree];
  for (int i = 0; i < indegree; ++i) {
    CALL(MPI_Irecv(static_cast<char*>(recvbuf) + rdispls[i] * recvwidth * width,
        (rdispls[i + 1] - rdispls[i]) * width, recvtype, sources[i], tag, comm,
        recvreqs + i));
  }
  for (int i = 0; i < outdegree; ++i) {
    CALL(MPI_Send(
        static_cast<char const*>(sendbuf) + sdispls[i] * sendwidth * width,
        (sdispls[i + 1] - sdispls[i]) * width, sendtype, destinations[i], tag,
        comm));
  }
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  delete[] recvreqs;
  end_code();
  return MPI_SUCCESS;
}

#endif  // end ifdef OMEGA_H_USE_MPI

template <typename T>
Read<T> Comm::allgather(T x) const {
#ifdef OMEGA_H_USE_MPI
  HostWrite<T> recvbuf(srcs_.size());
  CALL(Neighbor_allgather(host_srcs_, host_dsts_, &x, 1,
      MpiTraits<T>::datatype(), nonnull(recvbuf.data()), 1,
      MpiTraits<T>::datatype(), impl_));
  return recvbuf.write();
#else
  if (srcs_.size() == 1) return Read<T>({x});
  return Read<T>({});
#endif
}

template <typename T>
Read<T> Comm::alltoall(Read<T> x) const {
#ifdef OMEGA_H_USE_MPI
  HostWrite<T> recvbuf(srcs_.size());
  HostRead<T> sendbuf(x);
  CALL(Neighbor_alltoall(host_srcs_, host_dsts_, nonnull(sendbuf.data()), 1,
      MpiTraits<T>::datatype(), nonnull(recvbuf.data()), 1,
      MpiTraits<T>::datatype(), impl_));
  return recvbuf.write();
#else
  return x;
#endif
}

#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)

template <typename T>
Read<T> self_send_part1(LO self_dst, LO self_src, Read<T>* p_sendbuf,
    Read<LO>* p_sdispls, Read<LO>* p_rdispls, Int width, LO threshold) {
  Read<T> self_data;
  if (self_dst < 0) return self_data;
  OMEGA_H_CHECK(self_src >= 0);
  auto sendbuf = *p_sendbuf;
  auto sdispls = *p_sdispls;
  auto rdispls = *p_rdispls;
  auto begin = sdispls.get(self_dst) * width;
  auto end = sdispls.get(self_dst + 1) * width;
  auto self_size = end - begin;
  if (self_size == sendbuf.size()) {
    self_data = sendbuf;
    sendbuf = Read<T>({});
    sdispls = LOs({0, 0});
    auto recvcounts_w = deep_copy(get_degrees(rdispls));
    recvcounts_w.set(self_src, 0);
    auto recvcounts = LOs(recvcounts_w);
    rdispls = offset_scan(recvcounts);
  } else {
    if (self_size * sizeof(T) < size_t(threshold)) return self_data;
    auto self_data_w = Write<T>(self_size);
    auto other_data_w = Write<T>(sendbuf.size() - self_size);
    auto f = OMEGA_H_LAMBDA(LO i) {
      if (i < begin)
        other_data_w[i] = sendbuf[i];
      else if (i < end)
        self_data_w[i - begin] = sendbuf[i];
      else
        other_data_w[i - self_size] = sendbuf[i];
    };
    parallel_for(sendbuf.size(), f, "self_send_part1");
    self_data = self_data_w;
    sendbuf = other_data_w;
    auto sendcounts_w = deep_copy(get_degrees(sdispls));
    auto recvcounts_w = deep_copy(get_degrees(rdispls));
    sendcounts_w.set(self_dst, 0);
    recvcounts_w.set(self_src, 0);
    auto sendcounts = LOs(sendcounts_w);
    auto recvcounts = LOs(recvcounts_w);
    sdispls = offset_scan(sendcounts);
    rdispls = offset_scan(recvcounts);
  }
  *p_sendbuf = sendbuf;
  *p_sdispls = sdispls;
  *p_rdispls = rdispls;
  return self_data;
}

template <typename T>
void self_send_part2(Read<T> self_data, LO self_src, Read<T>* p_recvbuf,
    Read<LO> rdispls, Int width) {
  if (!self_data.exists()) return;
  auto recvbuf = *p_recvbuf;
  if (recvbuf.size() == 0) {
    recvbuf = self_data;
  } else {
    auto begin = rdispls.get(self_src) * width;
    auto self_size = self_data.size();
    auto end = begin + self_size;
    auto recvbuf_w = Write<T>(recvbuf.size() + self_size);
    auto f = OMEGA_H_LAMBDA(LO i) {
      if (i < begin)
        recvbuf_w[i] = recvbuf[i];
      else if (i < end)
        recvbuf_w[i] = self_data[i - begin];
      else
        recvbuf_w[i] = recvbuf[i - self_size];
    };
    parallel_for(recvbuf_w.size(), f, "self_send_part2");
    recvbuf = recvbuf_w;
  }
  *p_recvbuf = recvbuf;
}

#endif

template <typename T>
Read<T> Comm::alltoallv(Read<T> sendbuf_dev, Read<LO> sdispls_dev,
    Read<LO> rdispls_dev, Int width) const {
  begin_code("Comm::alltoallv");
#ifdef OMEGA_H_USE_MPI
#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  auto self_data = self_send_part1(self_dst_, self_src_, &sendbuf_dev,
      &sdispls_dev, &rdispls_dev, width, library_->self_send_threshold());
#endif
  HostRead<LO> sdispls(sdispls_dev);
  HostRead<LO> rdispls(rdispls_dev);
  OMEGA_H_CHECK(sendbuf_dev.size() == sdispls.last() * width);
  int nrecvd = rdispls.last() * width;
#if defined(OMEGA_H_USE_CUDA) && !defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  HostWrite<T> recvbuf(nrecvd);
  HostRead<T> sendbuf(sendbuf_dev);
  OMEGA_H_CHECK(recvbuf.size() == rdispls.last() * width);
  CALL(Neighbor_alltoallv(host_srcs_, host_dsts_, width,
      nonnull(sendbuf.data()), nonnull(sdispls.data()),
      MpiTraits<T>::datatype(), nonnull(recvbuf.data()),
      nonnull(rdispls.data()), MpiTraits<T>::datatype(), impl_));
  auto recvbuf_dev = Read<T>(recvbuf.write());
  self_send_part2(self_data, self_src_, &recvbuf_dev, rdispls_dev, width);
#else   // !defined(OMEGA_H_USE_CUDA) || defined(OMEGA_H_USE_CUDA_AWARE_MPI)
  Write<T> recvbuf_dev_w(nrecvd);
  OMEGA_H_CHECK(recvbuf_dev_w.size() == rdispls.last() * width);
  CALL(Neighbor_alltoallv(host_srcs_, host_dsts_, width,
      nonnull(sendbuf_dev.data()), nonnull(sdispls.data()),
      MpiTraits<T>::datatype(), nonnull(recvbuf_dev_w.data()),
      nonnull(rdispls.data()), MpiTraits<T>::datatype(), impl_));
  Read<T> recvbuf_dev = recvbuf_dev_w;
#endif  // !defined(OMEGA_H_USE_CUDA) || defined(OMEGA_H_USE_CUDA_AWARE_MPI)
#else   // !defined(OMEGA_H_USE_MPI)
  (void)sdispls_dev;
  (void)rdispls_dev;
  (void)width;
  auto recvbuf_dev = sendbuf_dev;
#endif  // !defined(OMEGA_H_USE_MPI)
  end_code();
  return recvbuf_dev;
}

void Comm::barrier() const {
  OMEGA_H_TIME_FUNCTION;
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Barrier(impl_));
#endif
}

#undef CALL

#define INST(T)                                                                \
  template T Comm::allreduce(T x, Omega_h_Op op) const;                        \
  template T Comm::exscan(T x, Omega_h_Op op) const;                           \
  template void Comm::bcast(T& x) const;                                       \
  template Read<T> Comm::allgather(T x) const;                                 \
  template Read<T> Comm::alltoall(Read<T> x) const;                            \
  template Read<T> Comm::alltoallv(                                            \
      Read<T> sendbuf, Read<LO> sdispls, Read<LO> rdispls, Int width) const;
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace Omega_h
