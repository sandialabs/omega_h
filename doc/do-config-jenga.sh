cmake /lore/dibanez/trilinos/omega_h2 \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_INSTALL_PREFIX=/lore/dibanez/trilinos/omega_h2/install \
-DKokkos_PREFIX=/lore/dibanez/trilinos/trilinos/install \
-DOSH_USE_KOKKOS=ON \
-DOSH_ONE_FILE=OFF \
-DBUILD_SHARED_LIBS=ON
