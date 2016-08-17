#ifndef COMM_HPP
#define COMM_HPP

#include "internal.hpp"

namespace osh {

#ifdef OMPI_MPI_H
/* OpenMPI defines MPI_UNWEIGHTED using (void*)
 * which causes compile errors with strict
 * compile options
 */
#define OMEGA_H_MPI_UNWEIGHTED reinterpret_cast<int*>(MPI_UNWEIGHTED)
#else
#define OMEGA_H_MPI_UNWEIGHTED MPI_UNWEIGHTED
#endif

#ifdef OMEGA_H_USE_MPI
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
#endif  // OMEGA_H_USE_MPI

static_assert(sizeof(int) == 4, "Comm assumes 32-bit int");

#ifdef OMEGA_H_USE_MPI
inline MPI_Op mpi_op(osh_op op) {
  switch (op) {
    case OMEGA_H_MIN:
      return MPI_MIN;
    case OMEGA_H_MAX:
      return MPI_MAX;
    case OMEGA_H_SUM:
      return MPI_SUM;
  };
  NORETURN(MPI_MIN);
}
#endif

}  // end namespace osh

#endif
