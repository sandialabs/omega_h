#!/bin/bash
cmake . \
-DCMAKE_CXX_COMPILER=$HOME/kokkos-install/bin/nvcc_wrapper \
-DCMAKE_INSTALL_PREFIX=$HOME/omega_h2-install \
-DKokkos_PREFIX=$HOME/kokkos-install \
-DOSH_USE_KOKKOS=True \
-DOSH_USE_CUDA=True
