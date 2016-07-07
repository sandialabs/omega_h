#!/bin/bash
cmake $HOME/src/omega_h2 \
-DCMAKE_CXX_COMPILER=$HOME/install/gcc/mpich-3.2/bin/mpicxx \
-DCMAKE_INSTALL_PREFIX=$HOME/install/gcc/omega_h2-kokkos \
-DKokkos_PREFIX=$HOME/install/gcc/Trilinos \
-DOSH_USE_KOKKOS=True \
-DBUILD_SHARED_LIBS=ON
