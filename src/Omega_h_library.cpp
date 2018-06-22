#include <Omega_h_config.h>
#include <Omega_h_control.hpp>
#include <Omega_h_stack.hpp>

#include <csignal>
#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#ifdef OMEGA_H_USE_DWARF
#define BACKWARD_HAS_DWARF 1
#endif
#include <backward/backward.hpp>
namespace backward {
cfile_streambuf::int_type cfile_streambuf::underflow() {
  return traits_type::eof();
}
}  // namespace backward

#include "Omega_h_cmdline.hpp"
#include "Omega_h_library.hpp"

namespace Omega_h {

char* max_memory_stacktrace = nullptr;

void print_stacktrace(std::ostream& out, int max_frames) {
  ::backward::StackTrace st;
  st.load_here(std::size_t(max_frames));
  ::backward::Printer p;
  p.print(st, out);
}

char const* Library::static_version() { return OMEGA_H_SEMVER; }

char const* Library::version() { return static_version(); }

char const* Library::static_commit_id() { return OMEGA_H_COMMIT; }

char const* Library::commit_id() { return static_commit_id(); }

char const* Library::static_configure_options() { return OMEGA_H_CMAKE_ARGS; }

char const* Library::configure_options() { return static_configure_options(); }

#if defined(__GNUC__) && (!defined(__clang__))
#define _GNU_SOURCE 1
#include <fenv.h>
static void enable_floating_point_exceptions() {
  feclearexcept(FE_ALL_EXCEPT);
  // FE_INEXACT inexact result: rounding was necessary to store the result of an
  // earlier floating-point operation sounds like the above would happen in
  // almost any floating point operation involving non-whole numbers ??? As for
  // underflow, there are plenty of cases where we will have things like ((a +
  // eps) - (a)) -> eps, where eps can be arbitrarily close to zero (usually it
  // would have been zero with infinite precision).
  feenableexcept(FE_ALL_EXCEPT - FE_INEXACT - FE_UNDERFLOW);
}
#elif defined(__x86_64__) || defined(_M_X64)
#include <xmmintrin.h>
// Intel system
static void enable_floating_point_exceptions() {
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
#ifdef __clang__
#pragma clang diagnostic pop
#endif
}
#else
static void enable_floating_point_exceptions() { Omega_h_fail("FPE enabled but not supported"); }
#endif

void Library::initialize(char const* head_desc, int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
    ,
    MPI_Comm comm_mpi
#endif
) {
  std::string lib_desc = OMEGA_H_SEMVER;
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
  OMEGA_H_CHECK(MPI_SUCCESS == MPI_Initialized(&mpi_is_init));
  if (!mpi_is_init) {
    OMEGA_H_CHECK(MPI_SUCCESS == MPI_Init(argc, argv));
    we_called_mpi_init = true;
  } else {
    we_called_mpi_init = false;
  }
  MPI_Comm world_dup;
  MPI_Comm_dup(comm_mpi, &world_dup);
  world_ = CommPtr(new Comm(this, world_dup));
#else
  world_ = CommPtr(new Comm(this, false, false));
  self_ = CommPtr(new Comm(this, false, false));
#endif
  Omega_h::CmdLine cmdline;
  cmdline.add_flag(
      "--osh-memory", "print amount and stacktrace of max memory use");
  cmdline.add_flag(
      "--osh-time", "print amount of time spend in certain functions");
  cmdline.add_flag("--osh-signal", "catch signals and print a stacktrace");
  cmdline.add_flag("--osh-fpe", "enable floating-point exceptions");
  cmdline.add_flag("--osh-silent", "suppress all output");
  auto& self_send_flag =
      cmdline.add_flag("--osh-self-send", "control self send threshold");
  self_send_flag.add_arg<int>("value");
  if (argc && argv) {
    OMEGA_H_CHECK(cmdline.parse(world_, argc, *argv));
  }
  if (cmdline.parsed("--osh-time")) {
    Omega_h::perf::global_singleton_history =
      new Omega_h::perf::History();
  }
  if (cmdline.parsed("--osh-fpe")) {
    enable_floating_point_exceptions();
  }
  bool should_protect = cmdline.parsed("--osh-signal");
  self_send_threshold_ = 1000 * 1000;
  if (cmdline.parsed("--osh-self-send")) {
    self_send_threshold_ = cmdline.get<int>("--osh-self-send", "value");
  }
  silent_ = cmdline.parsed("--osh-silent");
#ifdef OMEGA_H_USE_KOKKOSCORE
  if (!Kokkos::is_initialized()) {
    OMEGA_H_CHECK(argc != nullptr);
    OMEGA_H_CHECK(argv != nullptr);
    Kokkos::initialize(*argc, *argv);
    we_called_kokkos_init = true;
  } else {
    we_called_kokkos_init = false;
  }
#endif
  if (should_protect) Omega_h_protect();
}

Library::Library(Library const& other)
    : world_(other.world_),
      self_(other.self_)
#ifdef OMEGA_H_USE_MPI
      ,
      we_called_mpi_init(other.we_called_mpi_init)
#endif
#ifdef OMEGA_H_USE_KOKKOSCORE
      ,
      we_called_kokkos_init(other.we_called_kokkos_init)
#endif
{
}

Library::~Library() {
  if (Omega_h::perf::global_singleton_history) {
    if (world_->rank() == 0) {
      Omega_h::perf::print_top_down_and_bottom_up(
          *Omega_h::perf::global_singleton_history);
    }
    delete Omega_h::perf::global_singleton_history;
  }
#ifdef OMEGA_H_USE_KOKKOSCORE
  if (we_called_kokkos_init) {
    Kokkos::finalize();
    we_called_kokkos_init = false;
  }
#endif
  // need to destroy all Comm objects prior to MPI_Finalize()
  world_ = CommPtr();
  self_ = CommPtr();
#ifdef OMEGA_H_USE_MPI
  if (we_called_mpi_init) {
    OMEGA_H_CHECK(MPI_SUCCESS == MPI_Finalize());
    we_called_mpi_init = false;
  }
#endif
  delete[] Omega_h::max_memory_stacktrace;
}

CommPtr Library::world() { return world_; }

CommPtr Library::self() {
#ifdef OMEGA_H_USE_MPI
  if (!self_) {
    MPI_Comm self_dup;
    MPI_Comm_dup(MPI_COMM_SELF, &self_dup);
    self_ = CommPtr(new Comm(this, self_dup));
  }
#endif
  return self_;
}

LO Library::self_send_threshold() const { return self_send_threshold_; }

}  // end namespace Omega_h
