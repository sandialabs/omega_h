#include <cstdarg>
#include <cstdio>

#ifdef USE_KOKKOS
static bool we_called_kokkos_init = false;
#endif

void init(int& argc, char**& argv) {
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
}

void fail(char const* format, ...)
{
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  abort();
}
