#!/bin/bash
cmake /lore/dibanez/trilinos/src/omega_h2 \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_INSTALL_PREFIX=/lore/dibanez/trilinos/install/omega_h2 \
-DKokkos_PREFIX=/lore/dibanez/trilinos/install/Trilinos \
-DOmega_h_USE_Kokkos=ON \
-DOmega_h_ONE_FILE=OFF \
-DBUILD_SHARED_LIBS=ON
