#ifdef OMPI_MPI_H
/* OpenMPI defines MPI_UNWEIGHTED using (void*)
 * which causes compile errors with strict
 * compile options
 */
#define OSH_MPI_UNWEIGHTED reinterpret_cast<int*>(MPI_UNWEIGHTED)
#else
#define OSH_MPI_UNWEIGHTED MPI_UNWEIGHTED
#endif

#ifdef OSH_USE_MPI
template <class T>
struct MpiTraits;

template <>
struct MpiTraits<char> {
  static MPI_Datatype datatype() { return MPI_CHAR; }
};

template <>
struct MpiTraits<I8> {
  static MPI_Datatype datatype() { return MPI_INT8_T; }
};

template <>
struct MpiTraits<I32> {
  static MPI_Datatype datatype() { return MPI_INT32_T; }
};

template <>
struct MpiTraits<I64> {
  static MPI_Datatype datatype() { return MPI_INT64_T; }
};

template <>
struct MpiTraits<double> {
  static MPI_Datatype datatype() { return MPI_DOUBLE; }
};
#endif //OSH_USE_MPI

class Comm;

typedef std::shared_ptr<Comm> CommPtr;

enum ReduceOp {
  MIN,
  MAX,
  SUM
};

#ifdef OSH_USE_MPI
inline MPI_Op mpi_op(ReduceOp op) {
  switch (op) {
    case MIN: return MPI_MIN;
    case MAX: return MPI_MAX;
    case SUM: return MPI_SUM;
  };
  NORETURN(MPI_MIN);
}
#endif

static_assert(sizeof(int) == 4, "Comm assumes 32-bit int");

class Comm {
#ifdef OSH_USE_MPI
  MPI_Comm impl_;
#endif
  Read<I32> srcs_;
  HostRead<I32> host_srcs_;
  Read<I32> dsts_;
  HostRead<I32> host_dsts_;
public:
  Comm();
#ifdef OSH_USE_MPI
  Comm(MPI_Comm impl);
#else
  Comm(bool sends_to_self);
#endif
  ~Comm();
  static CommPtr world();
  static CommPtr self();
  I32 rank() const;
  I32 size() const;
  CommPtr dup() const;
  CommPtr split(I32 color, I32 key) const;
  CommPtr graph(Read<I32> dsts) const;
  CommPtr graph_adjacent(Read<I32> srcs, Read<I32> dsts) const;
  CommPtr graph_inverse() const;
  Read<I32> sources() const;
  Read<I32> destinations() const;
  template <typename T>
  T allreduce(T x, ReduceOp op) const;
  Int128 add_int128(Int128 x) const;
  template <typename T>
  T exscan(T x, ReduceOp op) const;
  template <typename T>
  void bcast(T& x) const;
  void bcast_string(std::string& s) const;
  template <typename T>
  Read<T> allgather(T x) const;
  template <typename T>
  Read<T> alltoall(Read<T> x) const;
  template <typename T>
  Read<T> alltoallv(Read<T> sendbuf,
      Read<LO> sendcounts, Read<LO> sdispls,
      Read<LO> recvcounts, Read<LO> rdispls) const;
  void barrier() const;
};
