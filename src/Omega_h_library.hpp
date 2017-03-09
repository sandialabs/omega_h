#ifndef OMEGA_H_LIBRARY_HPP
#define OMEGA_H_LIBRARY_HPP

#include <map>

#include <Omega_h_comm.hpp>

namespace Omega_h {

class Library {
 public:
  Library(Library const&);
  inline Library(int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
      ,
      MPI_Comm comm_mpi = MPI_COMM_WORLD
#endif
      ) {
    initialize(OMEGA_H_VERSION, argc, argv
#ifdef OMEGA_H_USE_MPI
        ,
        comm_mpi
#endif
        );
  }
  ~Library();
  static char const* static_version();
  char const* version();
  CommPtr world();
  CommPtr self();
  void add_to_timer(std::string const& name, double nsecs);
  LO self_send_threshold() const;
  bool should_time_;
  LO self_send_threshold_;
  bool silent_;

 private:
  void initialize(char const* head_desc, int* argc, char*** argv
#ifdef OMEGA_H_USE_MPI
      ,
      MPI_Comm comm_mpi
#endif
      );
  CommPtr world_;
  CommPtr self_;
#ifdef OMEGA_H_USE_MPI
  bool we_called_mpi_init;
#endif
#ifdef OMEGA_H_USE_KOKKOS
  bool we_called_kokkos_init;
#endif
  std::map<std::string, double> timers;
};

}

#endif
