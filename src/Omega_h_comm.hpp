#ifndef OMEGA_H_COMM_HPP
#define OMEGA_H_COMM_HPP

#include <memory>

#include <Omega_h_c.h>
#include <Omega_h_defines.hpp>
#include <Omega_h_int128.hpp>
#include <Omega_h_array.hpp>

namespace Omega_h {

static_assert(sizeof(int) == 4, "Omega_h::Comm assumes 32-bit int");

class Library;
class Comm;

typedef std::shared_ptr<Comm> CommPtr;

class Comm {
#ifdef OMEGA_H_USE_MPI
  MPI_Comm impl_;
#endif
  Library* library_;
  Read<I32> srcs_;
  Read<I32> dsts_;
  HostRead<I32> host_srcs_;
  HostRead<I32> host_dsts_;
  LO self_src_;
  LO self_dst_;

 public:
  Comm();
#ifdef OMEGA_H_USE_MPI
  Comm(Library* library, MPI_Comm impl);
#else
  Comm(Library* library, bool is_graph, bool sends_to_self);
#endif
  ~Comm();
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
  T allreduce(T x, Omega_h_Op op) const;
  bool reduce_or(bool x) const;
  bool reduce_and(bool x) const;
  Int128 add_int128(Int128 x) const;
  template <typename T>
  T exscan(T x, Omega_h_Op op) const;
  template <typename T>
  void bcast(T& x) const;
  void bcast_string(std::string& s) const;
  template <typename T>
  Read<T> allgather(T x) const;
  template <typename T>
  Read<T> alltoall(Read<T> x) const;
  template <typename T>
  Read<T> alltoallv(Read<T> sendbuf, Read<LO> sendcounts, Read<LO> sdispls,
      Read<LO> recvcounts, Read<LO> rdispls) const;
  void barrier() const;
};

#ifdef OMEGA_H_USE_MPI

#ifdef OMPI_MPI_H
/* OpenMPI defines MPI_UNWEIGHTED using (void*)
 * which causes compile errors with strict
 * compile options
 */
#define OMEGA_H_MPI_UNWEIGHTED reinterpret_cast<int*>(MPI_UNWEIGHTED)
#else
#define OMEGA_H_MPI_UNWEIGHTED MPI_UNWEIGHTED
#endif

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

template <>
struct MpiTraits<unsigned> {
  static MPI_Datatype datatype() { return MPI_UNSIGNED; }
};

template <>
struct MpiTraits<unsigned long> {
  static MPI_Datatype datatype() { return MPI_UNSIGNED_LONG; }
};

inline MPI_Op mpi_op(Omega_h_Op op) {
  switch (op) {
    case OMEGA_H_MIN:
      return MPI_MIN;
    case OMEGA_H_MAX:
      return MPI_MAX;
    case OMEGA_H_SUM:
      return MPI_SUM;
  };
  OMEGA_H_NORETURN(MPI_MIN);
}
#endif

}  // end namespace Omega_h

#endif
