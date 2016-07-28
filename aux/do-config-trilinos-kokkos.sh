#!/bin/bash
cmake $HOME/src/Trilinos \
 -Wno-dev \
 -DCMAKE_INSTALL_PREFIX:PATH=$HOME/install/gcc/Trilinos-kokkos \
 -DCMAKE_BUILD_TYPE:STRING=NONE \
 -DBUILD_SHARED_LIBS:BOOL=ON \
 -DTPL_ENABLE_MPI:BOOL=ON \
 -DMPI_BASE_DIR:PATH=$HOME/install/gcc/mpich-3.2 \
 -DCMAKE_CXX_COMPILER:FILEPATH=$HOME/install/gcc/mpich-3.2/bin/mpicxx \
 -DCMAKE_C_COMPILER:FILEPATH=$HOME/install/gcc/mpich-3.2/bin/mpicc \
 -DTrilinos_ENABLE_Fortran=OFF \
 -DCMAKE_CXX_FLAGS:STRING='-O2 -g' \
 -DCMAKE_C_FLAGS:STRING='-O2 -g' \
 -DTrilinos_CXX11_FLAGS:STRING='-std=c++11 -Wno-unused-local-typedefs -Wno-sign-compare' \
 -DTrilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
 -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
 -DTrilinos_ENABLE_Kokkos:BOOL=ON \
 -DTrilinos_ENABLE_KokkosCore:BOOL=ON \
 -DTrilinos_ENABLE_KokkosContainers:BOOL=ON \
 -DTrilinos_ENABLE_KokkosExample:BOOL=ON \
2>&1 | tee config_log
