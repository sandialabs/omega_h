#!/bin/bash -ex
./unit_tests
mpirun -np 4 ./mpi_tests
mpirun -np 2 ./manual_test
