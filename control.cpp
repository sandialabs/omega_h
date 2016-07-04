#ifdef OSH_USE_MPI
static bool we_called_mpi_init = false;
#endif
#ifdef OSH_USE_KOKKOS
static bool we_called_kokkos_init = false;
#endif

extern "C" void osh_init_internal(int* argc, char*** argv,
                                  char const* head_desc) {
  std::string lib_desc = OSH_VERSION;
  if (lib_desc != head_desc) {
    std::stringstream msg;
    msg << "omega_h description string mismatch.\n";
    msg << "header says: " << head_desc << '\n';
    msg << "library says: " << lib_desc << '\n';
    std::string msg_str = msg.str();
    osh_fail("%s\n", msg_str.c_str());
  }
#ifdef OSH_USE_MPI
  int mpi_is_init;
  CHECK(MPI_SUCCESS == MPI_Initialized(&mpi_is_init));
  if (!mpi_is_init) {
    CHECK(MPI_SUCCESS == MPI_Init(argc, argv));
    we_called_mpi_init = true;
  }
#endif
#ifdef OSH_USE_KOKKOS
  if (!Kokkos::DefaultExecutionSpace::is_initialized()) {
    Kokkos::initialize(*argc, *argv);
    we_called_kokkos_init = true;
  }
#endif
  (void)argc;
  (void)argv;
  protect();
}

extern "C" void osh_finalize(void) {
#ifdef OSH_USE_KOKKOS
  if (we_called_kokkos_init) {
    Kokkos::finalize();
    we_called_kokkos_init = false;
  }
#endif
#ifdef OSH_USE_MPI
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

extern "C" void osh_fail(char const* format, ...) {
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  abort();
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

Library::~Library() { osh_finalize(); }

CommPtr Library::world() const { return Comm::world(); }

CommPtr Library::self() const { return Comm::self(); }
