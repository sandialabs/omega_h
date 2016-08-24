#!/bin/bash
MPI_BASE_DIR=/opt/intel/compilers_and_libraries/linux/mpi/intel64
ZLIB_DIR=/usr/local/mic
cmake $HOME/src/omega_h2 \
-DCMAKE_CXX_COMPILER:FILEPATH=${MPI_BASE_DIR}/bin/mpiicpc \
-DCMAKE_INSTALL_PREFIX:PATH=$HOME/install/omega_h \
-DZLIB_ROOT:PATH=${ZLIB_DIR} \
-DOmega_h_USE_Kokkos:BOOL=True \
-DKokkos_PREFIX:PATH=$HOME/install/Trilinos \
-DOmega_h_USE_OpenMP:BOOL=True \
-DOmega_h_EXTRA_CXX_FLAGS:STRING="-w -mmic -DPREC_TIMER -restrict -fasm-blocks -DDEVICE=1wq" \
-DCMAKE_VERBOSE_MAKEFILE=OFF \
-DOmega_h_ONE_FILE=OFF \
2>&1 | tee config_log
