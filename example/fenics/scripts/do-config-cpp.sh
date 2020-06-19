#!/bin/bash -ex
cmake $HOME/shared/omega_h-dolfin \
-DCMAKE_CXX_COMPILER:FILEPATH=`which mpicxx` \
-DCMAKE_C_COMPILER:FILEPATH=`which mpicc` \
-DCMAKE_BUILD_TYPE:STRING="" \
-DCMAKE_CXX_FLAGS:STRING="-g -O3 -std=c++11 -fno-omit-frame-pointer" \
-DOMEGA_H_PREFIX="$HOME/shared/omega_h/install" \
2>&1 | tee config_log
