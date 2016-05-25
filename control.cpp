#include <cstdarg>
#include <cstdio>

#ifdef USE_MPI
static bool we_called_mpi_init = false;
#endif
#ifdef USE_KOKKOS
static bool we_called_kokkos_init = false;
#endif

void init(int& argc, char**& argv) {
#ifdef USE_MPI
  int mpi_is_init;
  CHECK(MPI_SUCCESS == MPI_Initialized(&mpi_is_init));
  if (!mpi_is_init) {
    CHECK(MPI_SUCCESS == MPI_Init(&argc, &argv));
    we_called_mpi_init = true;
  }
#endif
#ifdef USE_KOKKOS
  if (!Kokkos::DefaultExecutionSpace::is_initialized()) {
    Kokkos::initialize(argc, argv);
    we_called_kokkos_init = true;
  }
#endif
  (void)argc;
  (void)argv;
}

void fini() {
#ifdef USE_KOKKOS
  if (we_called_kokkos_init) {
    Kokkos::finalize();
    we_called_kokkos_init = false;
  }
#endif
#ifdef USE_MPI
  if (we_called_mpi_init) {
    CHECK(MPI_SUCCESS == MPI_Finalize());
  }
#endif
}

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-nonliteral"
#endif

void fail(char const* format, ...)
{
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  abort();
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
