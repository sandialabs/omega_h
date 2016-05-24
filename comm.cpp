#define CALL(f) CHECK(MPI_SUCCESS == (f))

Comm::Comm() {
#ifdef USE_MPI
  impl_ = MPI_COMM_NULL;
#endif
}

#ifdef USE_MPI
Comm::Comm(MPI_Comm impl):impl_(impl) {
  int topo_type;
  CALL(MPI_Topo_test(impl, &topo_type));
  if (topo_type == MPI_DIST_GRAPH) {
    int nin, nout, is_weighted;
    CALL(MPI_Dist_graph_neighbors_count(impl, &nin, &nout, &is_weighted));
    HostWrite<I32> sources(nin);
    HostWrite<I32> destinations(nin);
    CALL(MPI_Dist_graph_neighbors(
          impl,
          nin,
          &sources[0],
          OSH_MPI_UNWEIGHTED,
          nout,
          &destinations[0],
          OSH_MPI_UNWEIGHTED));
    srcs_ = sources.write();
    dsts_ = destinations.write();
    host_srcs_ = HostRead<I32>(srcs);
    host_dsts_ = HostRead<I32>(dsts);
  }
}
#endif

Comm::~Comm() {
  CALL(MPI_Comm_free(&impl_));
}

static CommPtr Comm::world() {
  MPI_Comm impl;
  CALL(MPI_Comm_dup(MPI_COMM_WORLD, &impl));
  return new Comm(impl);
}

static CommPtr Comm::self() {
  MPI_Comm impl;
  CALL(MPI_Comm_dup(MPI_COMM_SELF, &impl));
  return new Comm(impl);
}

I32 Comm::rank() const {
  I32 r;
  CALL(MPI_Comm_rank(impl_, &r));
  return r;
}

I32 Comm::size() const {
  I32 s;
  CALL(MPI_Comm_size(impl_, &s));
  return s;
}

CommPtr Comm::dup() const {
  MPI_Comm impl2;
  CALL(MPI_Comm_dup(impl_, &impl2));
  return new Comm(impl2);
}

CommPtr Comm::split(I32 color, I32 key) const {
  MPI_Comm impl2;
  CALL(MPI_Comm_split(impl_, color, key, &impl2));
  return impl2;
}

CommPtr Comm::graph(Read<I32> dsts) const {
  MPI_Comm impl2;
  int n = 1;
  int sources[1] = {rank()};
  int degrees[1] = {dsts.size()};
  HostRead<I32> destinations(dsts);
  int reorder = 0;
  CALL(MPI_Dist_graph_create(
        impl_,
        n,
        sources,
        degrees,
        &destinations[0],
        OSH_MPI_UNWEIGHTED,
        MPI_INFO_NULL,
        reorder,
        &impl2));
  return new Comm(impl2);
}

CommPtr Comm::graph_adjacent(Read<I32> srcs, Read<I32> dsts) const {
  MPI_Comm impl2;
  HostRead<I32> sources(srcs);
  HostRead<I32> destinations(dsts);
  int reorder = 0;
  CALL(MPI_Dist_graph_create_adjacent(
        impl_,
        sources.size(), &sources[0],
        OSH_MPI_UNWEIGHTED,
        destinations.size(), &destinations[0],
        OSH_MPI_UNWEIGHTED,
        MPI_INFO_NULL, reorder, &impl2));
  return new Comm(impl2);
}

Read<I32> Comm::sources() const {
  return srcs_;
}

Read<I32> Comm::destinations() const {
  return dsts_;
}

template <typename T>
T Comm::allreduce(T x, ReduceOp op) const {
  CALL(MPI_Allreduce(MPI_IN_PLACE, &x, 1, MpiTraits<T>::datatype(),
        mpi_op(op), impl_));
  return x;
}

template <typename T>
T Comm::exscan(T x, ReduceOp op) const {
  CALL(MPI_Exscan(MPI_IN_PLACE, &x, 1, MpiTraits<T>::datatype(),
        mpi_op(op), impl_));
  if (!rank())
    return 0;
  return value;
}

template <typename T>
void Comm::bcast(T& x) const {
  CALL(MPI_Bcast(&x, 1, MpiTraits<T>::datatype(), 0, impl_));
}

void Comm::bcast_string(std::string& s) const {
  I32 len = static_cast<I32>(s.length());
  len = bcast(len);
  s.resize(len);
  CALL(MPI_Bcast(s.c_str(), len, MPI_CHAR, 0, impl_));
}

/* custom implementation of MPI_Neighbor_allgather
 * in the case that we are using an MPI older
 * than version 3.0
 */

static int Neighbor_allgather(
    HostRead<I32> sources,
    HostRead<I32> destinations,
    const void *sendbuf,
    int sendcount,
    MPI_Datatype sendtype,
    void *recvbuf,
    int recvcount,
    MPI_Datatype recvtype,
    MPI_Comm comm)
{
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
    CALL(MPI_Irecv(
          static_cast<char*>(recvbuf) + i * recvwidth,
          recvcount, recvtype, sources[i], tag, comm,
          recvreqs + i));
  delete [] sources;
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(
          sendbuf,
          sendcount, sendtype, destinations[i], tag, comm,
          sendreqs + i));
  delete [] destinations;
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  delete [] sendreqs;
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  delete [] recvreqs;
  return MPI_SUCCESS;
#else
  return MPI_Neighbor_allgather(
      sendbuf, sendcount, sendtype,
      recvbuf, recvcount, recvtype,
      comm);
#endif
}

/* custom implementation of MPI_Neighbor_alltoall
 * in the case that we are using an MPI older
 * than version 3.0
 */

static int Neighbor_alltoall(
    HostRead<I32> sources,
    HostRead<I32> destinations,
    const void *sendbuf,
    int sendcount,
    MPI_Datatype sendtype,
    void *recvbuf,
    int recvcount,
    MPI_Datatype recvtype,
    MPI_Comm comm)
{
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
    CALL(MPI_Irecv(
          static_cast<char*>(recvbuf) + i * recvwidth,
          recvcount, recvtype, sources[i], tag, comm,
          recvreqs + i));
  delete [] sources;
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(
          static_cast<char const*>(sendbuf) + i * sendwidth,
          sendcount, sendtype, destinations[i], tag, comm,
          sendreqs + i));
  delete [] destinations;
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  delete [] sendreqs;
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  delete [] recvreqs;
  return MPI_SUCCESS;
#else
  return MPI_Neighbor_alltoall(
      sendbuf, sendcount, sendtype,
      recvbuf, recvcount, recvtype,
      comm);
#endif
}

/* custom implementation of MPI_Neighbor_alltoallv
 * in the case that we are using an MPI older
 * than version 3.0
 */

static int Neighbor_alltoallv(
    HostRead<I32> sources,
    HostRead<I32> destinations,
    const void *sendbuf,
    const int sendcounts[],
    const int sdispls[],
    MPI_Datatype sendtype,
    void *recvbuf,
    const int recvcounts[],
    const int rdispls[],
    MPI_Datatype recvtype,
    MPI_Comm comm)
{
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
    CALL(MPI_Irecv(
          static_cast<char*>(recvbuf) + rdispls[i] * recvwidth,
          recvcounts[i], recvtype, sources[i], tag, comm,
          recvreqs + i));
  delete [] sources;
  CALL(MPI_Barrier(comm));
  for (int i = 0; i < outdegree; ++i)
    CALL(MPI_Isend(
          static_cast<char const*>(sendbuf) + sdispls[i] * sendwidth,
          sendcounts[i], sendtype, destinations[i], tag, comm,
          sendreqs + i));
  delete [] destinations;
  CALL(MPI_Waitall(outdegree, sendreqs, MPI_STATUSES_IGNORE));
  delete [] sendreqs;
  CALL(MPI_Waitall(indegree, recvreqs, MPI_STATUSES_IGNORE));
  delete [] recvreqs;
  return MPI_SUCCESS;
#else
  return MPI_Neighbor_alltoallv(
      sendbuf, sendcounts, sdispls, sendtype,
      recvbuf, recvcounts, rdispls, recvtype,
      comm);
#endif
}

template <typename T>
Read<T> Comm::allgather(T x) const {
  HostWrite<T> recvbuf(srcs_.size());
  CALL(Neighbor_allgather(
        host_srcs_, host_dsts_,
        &x,          1, MpiTraits<T>::datatype(),
        &recvbuf[0], 1, MpiTraits<T>::datatype(),
        impl_));
  return recvbuf.write();
}

template <typename T>
Read<T> Comm::alltoall(Read<T> x) const {
  HostWrite<T> recvbuf(srcs_.size());
  HostRead<T,int> sendbuf(x);
  CALL(Neighbor_alltoall(
        host_srcs_, host_dsts_,
        &sendbuf[0], 1, MpiTraits<T>::datatype(),
        &recvbuf[0], 1, MpiTraits<T>::datatype(),
        impl_));
  return recvbuf.write();
}

template <typename T>
Read<T> Comm::alltoallv(Read<T> sendbuf_dev,
    Read<LO> sendcounts_dev, Read<LO> sdispls_dev,
    Read<LO> recvcounts_dev, Read<LO> rdispls_dev) const {
  HostRead<T> sendbuf(sendbuf_dev);
  HostRead<LO> sendcounts(sendcounts_dev);
  HostRead<LO> recvcounts(recvcounts_dev);
  HostRead<LO> sdispls(sdispls_dev);
  HostRead<LO> rdispls(rdispls_dev);
  CHECK(rdispls.size() == recvcounts.size() + 1);
  int nrecvd = rdispls[recvcounts.size()];
  HostWrite<T> recvbuf(nrecvd);
  CALL(Neighbor_alltoallv(
        host_srcs_, host_dsts_,
        &sendbuf[0], &sendcounts[0], &sdispls[0], MpiTraits<T>::datatype(),
        &recvbuf[0], &recvcounts[0], &rdispls[0], MpiTraits<T>::datatype(),
        impl_));
  return recvbuf.write();
}

#undef CALL

#define INST_T(T) \
template T Comm::allreduce(T x, ReduceOp op) const; \
template T Comm::exscan(T x, ReduceOp op) const; \
template T Comm::bcast(T x) const; \
template Read<T> Comm::allgather(T x) const; \
template Read<T> Comm::alltoall(Read<T> x) const; \
template Read<T> Comm::alltoallv(Read<T> sendbuf, \
    Read<LO> sendcounts, Read<LO> sdispls, \
    Read<LO> recvcounts, Read<LO> rdispls) const;
INST_T(I8)
INST_T(I32)
INST_T(I64)
INST_T(Real)
#undef INST_T
