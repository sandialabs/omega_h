#include "control.hpp"

#include <cstdarg>
#include <sstream>
#include <iostream>

#include "comm.hpp"
#include "internal.hpp"
#include "protect.hpp"

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

namespace Omega_h {

bool should_log_memory = false;
char* max_memory_stacktrace = nullptr;

static bool remove_flag(int* argc, char*** argv, std::string const& flag) {
  for (int i = 0; i < *argc; ++i) {
    if (flag == (*argv)[i]) {
      --(*argc);
      for (int j = i; j < *argc; ++j) (*argv)[j] = (*argv)[j + 1];
      return true;
    }
  }
  return false;
}

void Library::initialize(char const* head_desc, int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
      , MPI_Comm comm_mpi
#endif
      ) {
  std::string lib_desc = OMEGA_H_VERSION;
  if (lib_desc != head_desc) {
    std::stringstream msg;
    msg << "omega_h description string mismatch.\n";
    msg << "header says: " << head_desc << '\n';
    msg << "library says: " << lib_desc << '\n';
    std::string msg_str = msg.str();
    Omega_h_fail("%s\n", msg_str.c_str());
  }
  Omega_h::should_log_memory = remove_flag(argc, argv, "--osh-log-mem");
#ifdef OMEGA_H_USE_MPI
  int mpi_is_init;
  CHECK(MPI_SUCCESS == MPI_Initialized(&mpi_is_init));
  if (!mpi_is_init) {
    CHECK(MPI_SUCCESS == MPI_Init(argc, argv));
    we_called_mpi_init = true;
  }
  MPI_Comm world_dup;
  MPI_Comm_dup(MPI_COMM_WORLD, &world_dup);
  world_ = CommPtr(new Comm(world_dup));
#else
  world_ = CommPtr(new Comm());
  self_ = CommPtr(new Comm());
#endif
#ifdef OMEGA_H_USE_KOKKOS
  if (!Kokkos::DefaultExecutionSpace::is_initialized()) {
    CHECK(argc != nullptr);
    CHECK(argv != nullptr);
    Kokkos::initialize(*argc, *argv);
    we_called_kokkos_init = true;
  }
#endif
#ifdef OMEGA_H_PROTECT
  protect();
#endif
}

Library::Library(Library const& other):
  world_(other.world_),
  self_(other.self_)
#ifdef OMEGA_H_USE_MPI
  ,we_called_mpi_init(other.we_called_mpi_init)
#endif
#ifdef OMEGA_H_USE_KOKKOS
  ,we_called_kokkos_init(other.we_called_kokkos_init)
#endif
{
}

Library::~Library() {
#ifdef OMEGA_H_USE_KOKKOS
  if (we_called_kokkos_init) {
    Kokkos::finalize();
    we_called_kokkos_init = false;
  }
#endif
  if (Omega_h::should_log_memory) {
    auto mem_used = get_max_bytes();
    auto max_mem_used = std::size_t(
        world_->allreduce(I64(mem_used), OMEGA_H_MAX));
    auto max_mem_rank = (mem_used == max_mem_used) ?
      world_->rank() : world_->size();
    max_mem_rank = world_->allreduce(max_mem_rank, OMEGA_H_MIN);
    if (world_->rank() == max_mem_rank) {
      std::cout << "maximum Omega_h memory usage: " << mem_used << '\n';
      std::cout << Omega_h::max_memory_stacktrace;
    }
  }
#ifdef OMEGA_H_USE_MPI
  if (we_called_mpi_init) {
    CHECK(MPI_SUCCESS == MPI_Finalize());
    we_called_mpi_init = false;
  }
#endif
  delete [] Omega_h::max_memory_stacktrace;
}

CommPtr Library::world() const { return world_; }

CommPtr Library::self() const {
#ifdef OMEGA_H_USE_MPI
  if (!self_) {
    MPI_Comm self_dup;
    MPI_Comm_dup(MPI_COMM_WORLD, &self_dup);
    self_ = CommPtr(new Comm(self_dup));
  }
#endif
  return self_;
}

}  // end namespace Omega_h
