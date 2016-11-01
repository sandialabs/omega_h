#include "comm.hpp"
#include "internal.hpp"
#include "protect.hpp"

#include <cstdarg>
#include <sstream>

namespace Omega_h {

#ifdef OMEGA_H_USE_MPI
static bool we_called_mpi_init = false;
#endif
#ifdef OMEGA_H_USE_KOKKOS
static bool we_called_kokkos_init = false;
#endif

extern "C" void Omega_h_init_internal(
    int* argc, char*** argv, char const* head_desc) {
  std::string lib_desc = OMEGA_H_VERSION;
  if (lib_desc != head_desc) {
    std::stringstream msg;
    msg << "omega_h description string mismatch.\n";
    msg << "header says: " << head_desc << '\n';
    msg << "library says: " << lib_desc << '\n';
    std::string msg_str = msg.str();
    Omega_h_fail("%s\n", msg_str.c_str());
  }
#ifdef OMEGA_H_USE_MPI
  int mpi_is_init;
  CHECK(MPI_SUCCESS == MPI_Initialized(&mpi_is_init));
  if (!mpi_is_init) {
    CHECK(MPI_SUCCESS == MPI_Init(argc, argv));
    we_called_mpi_init = true;
  }
#endif
#ifdef OMEGA_H_USE_KOKKOS
  if (!Kokkos::DefaultExecutionSpace::is_initialized()) {
    CHECK(argc != nullptr);
    CHECK(argv != nullptr);
    Kokkos::initialize(*argc, *argv);
    we_called_kokkos_init = true;
  }
#endif
  (void)argc;
  (void)argv;
#ifdef OMEGA_H_PROTECT
  protect();
#endif
}

extern "C" void Omega_h_finalize(void) {
#ifdef OMEGA_H_USE_KOKKOS
  if (we_called_kokkos_init) {
    Kokkos::finalize();
    we_called_kokkos_init = false;
  }
#endif
#ifndef OMEGA_H_USE_KOKKOS
  std::size_t mem_used = Omega_h::get_max_bytes();
  std::size_t max_mem_used = mem_used;
#ifdef OMEGA_H_USE_MPI
  MPI_Reduce(&mem_used, &max_mem_used, 1,
      Omega_h::MpiTraits<std::size_t>::datatype(), MPI_MAX, 0, MPI_COMM_WORLD);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
#endif
    printf("maximum Omega_h memory usage: %zu bytes\n", mem_used);
#ifdef OMEGA_H_USE_MPI
  }
#endif
#endif
#ifdef OMEGA_H_USE_MPI
  if (we_called_mpi_init) {
    CHECK(MPI_SUCCESS == MPI_Finalize());
    we_called_mpi_init = false;
  }
#endif
}

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

extern "C" void Omega_h_fail(char const* format, ...) {
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  abort();
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

Library::~Library() { Omega_h_finalize(); }

CommPtr Library::world() const { return Comm::world(); }

CommPtr Library::self() const { return Comm::self(); }

}  // end namespace Omega_h
