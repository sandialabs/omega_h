#include "comm.hpp"

#include "array.hpp"
#include "int128.hpp"

namespace osh {

#define CALL(f) CHECK(MPI_SUCCESS == (f))

Comm::Comm() {
#ifdef OMEGA_H_USE_MPI
  impl_ = MPI_COMM_NULL;
#endif
}

#ifdef OMEGA_H_USE_MPI
Comm::Comm(MPI_Comm impl) : impl_(impl) {
  int topo_type;
  CALL(MPI_Topo_test(impl, &topo_type));
  if (topo_type == MPI_DIST_GRAPH) {
    int nin, nout, is_weighted;
    CALL(MPI_Dist_graph_neighbors_count(impl, &nin, &nout, &is_weighted));
    HostWrite<I32> sources(nin);
    HostWrite<I32> destinations(nout);
    CALL(MPI_Dist_graph_neighbors(impl, nin, sources.data(), OMEGA_H_MPI_UNWEIGHTED,
        nout, destinations.data(), OMEGA_H_MPI_UNWEIGHTED));
    srcs_ = sources.write();
    dsts_ = destinations.write();
    host_srcs_ = HostRead<I32>(srcs_);
    host_dsts_ = HostRead<I32>(dsts_);
  }
}
#else
Comm::Comm(bool is_graph, bool sends_to_self) {
  if (is_graph) {
    if (sends_to_self) {
      srcs_ = Read<LO>({0});
    } else {
      srcs_ = Read<LO>({});
    }
    dsts_ = srcs_;
    host_srcs_ = HostRead<I32>(srcs_);
    host_dsts_ = HostRead<I32>(dsts_);
  } else {
    CHECK(!sends_to_self);
  }
}
#endif

Comm::~Comm() {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Comm_free(&impl_));
#endif
}

CommPtr Comm::world() {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl;
  CALL(MPI_Comm_dup(MPI_COMM_WORLD, &impl));
  return CommPtr(new Comm(impl));
#else
  return CommPtr(new Comm());
#endif
}

CommPtr Comm::self() {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl;
  CALL(MPI_Comm_dup(MPI_COMM_SELF, &impl));
  return CommPtr(new Comm(impl));
#else
  return CommPtr(new Comm());
#endif
}

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
  return CommPtr(new Comm(impl2));
#else
  return CommPtr(new Comm(srcs_.exists(), srcs_.exists() && srcs_.size() == 1));
#endif
}

CommPtr Comm::split(I32 color, I32 key) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  CALL(MPI_Comm_split(impl_, color, key, &impl2));
  return CommPtr(new Comm(impl2));
#else
  (void)color;
  (void)key;
  return CommPtr(new Comm());
#endif
}

CommPtr Comm::graph(Read<I32> dsts) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  int n = 1;
  int sources[1] = {rank()};
  int degrees[1] = {dsts.size()};
  HostRead<I32> destinations(dsts);
  int reorder = 0;
  CALL(MPI_Dist_graph_create(impl_, n, sources, degrees, destinations.data(),
      OMEGA_H_MPI_UNWEIGHTED, MPI_INFO_NULL, reorder, &impl2));
  return CommPtr(new Comm(impl2));
#else
  return CommPtr(new Comm(true, dsts.size() == 1));
#endif
}

CommPtr Comm::graph_adjacent(Read<I32> srcs, Read<I32> dsts) const {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl2;
  HostRead<I32> sources(srcs);
  HostRead<I32> destinations(dsts);
  int reorder = 0;
  CALL(MPI_Dist_graph_create_adjacent(impl_, sources.size(), sources.data(),
      OMEGA_H_MPI_UNWEIGHTED, destinations.size(), destinations.data(),
      OMEGA_H_MPI_UNWEIGHTED, MPI_INFO_NULL, reorder, &impl2));
  return CommPtr(new Comm(impl2));
#else
  CHECK(srcs == dsts);
  return CommPtr(new Comm(true, dsts.size() == 1));
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
 * in the case that we are using an MPI older
 * than version 3.0
 */

static int Neighbor_alltoallv(HostRead<I32> sources, HostRead<I32> destinations,
    const void* sendbuf, const int sendcounts[], const int sdispls[],
    MPI_Datatype sendtype, void* recvbuf, const int recvcounts[],
    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm) {
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
    CALL(MPI_Irecv(static_cast<char*>(recvbuf) + rdispls[i] * recvwidth,
        recvcounts[i], recvtype, sources[i], tag, comm, recvreqs + i));
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(static_cast<char const*>(sendbuf) + sdispls[i] * sendwidth,
        sendcounts[i], sendtype, destinations[i], tag, comm, sendreqs + i));
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  delete[] sendreqs;
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  delete[] recvreqs;
  return MPI_SUCCESS;
#else
  (void)sources;
  (void)destinations;
  return MPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf,
      recvcounts, rdispls, recvtype, comm);
#endif  // end if MPI_VERSION < 3
}

#endif  // end ifdef OMEGA_H_USE_MPI

template <typename T>
Read<T> Comm::allgather(T x) const {
#ifdef OMEGA_H_USE_MPI
  HostWrite<T> recvbuf(srcs_.size());
  CALL(Neighbor_allgather(host_srcs_, host_dsts_, &x, 1,
      MpiTraits<T>::datatype(), recvbuf.data(), 1, MpiTraits<T>::datatype(),
      impl_));
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
  CALL(Neighbor_alltoall(host_srcs_, host_dsts_, sendbuf.data(), 1,
      MpiTraits<T>::datatype(), recvbuf.data(), 1, MpiTraits<T>::datatype(),
      impl_));
  return recvbuf.write();
#else
  return x;
#endif
}

template <typename T>
Read<T> Comm::alltoallv(Read<T> sendbuf_dev, Read<LO> sendcounts_dev,
    Read<LO> sdispls_dev, Read<LO> recvcounts_dev, Read<LO> rdispls_dev) const {
#ifdef OMEGA_H_USE_MPI
  HostRead<T> sendbuf(sendbuf_dev);
  HostRead<LO> sendcounts(sendcounts_dev);
  HostRead<LO> recvcounts(recvcounts_dev);
  HostRead<LO> sdispls(sdispls_dev);
  HostRead<LO> rdispls(rdispls_dev);
  CHECK(rdispls.size() == recvcounts.size() + 1);
  int nrecvd = rdispls.last();
  HostWrite<T> recvbuf(nrecvd);
  CHECK(sendcounts.size() == host_dsts_.size());
  CHECK(recvcounts.size() == host_srcs_.size());
  CHECK(sdispls.size() == sendcounts.size() + 1);
  CHECK(sendbuf.size() == sdispls.last());
  CALL(Neighbor_alltoallv(host_srcs_, host_dsts_, sendbuf.data(),
      sendcounts.data(), sdispls.data(), MpiTraits<T>::datatype(),
      recvbuf.data(), recvcounts.data(), rdispls.data(),
      MpiTraits<T>::datatype(), impl_));
  return recvbuf.write();
#else
  (void)sendcounts_dev;
  (void)recvcounts_dev;
  (void)sdispls_dev;
  (void)rdispls_dev;
  return sendbuf_dev;
#endif
}

void Comm::barrier() const {
#ifdef OMEGA_H_USE_MPI
  CALL(MPI_Barrier(impl_));
#endif
}

#undef CALL

#define INST(T)                                                                \
  template T Comm::allreduce(T x, Omega_h_Op op) const;                            \
  template T Comm::exscan(T x, Omega_h_Op op) const;                               \
  template void Comm::bcast(T& x) const;                                       \
  template Read<T> Comm::allgather(T x) const;                                 \
  template Read<T> Comm::alltoall(Read<T> x) const;                            \
  template Read<T> Comm::alltoallv(Read<T> sendbuf, Read<LO> sendcounts,       \
      Read<LO> sdispls, Read<LO> recvcounts, Read<LO> rdispls) const;
INST(I8)
INST(I32)
INST(I64)
INST(Real)
#undef INST

}  // end namespace osh
