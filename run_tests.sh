#!/bin/bash -ex
mpirun -np 1 ./unit_tests
mpirun -np 4 ./mpi_tests
mpirun -np 1 ./ring_test
mpirun -np 1 ./warp_test
mpirun -np 2 ./warp_test
