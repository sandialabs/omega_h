#!/bin/bash -ex
mpirun -np 1 ./unit_tests
mpirun -np 4 ./mpi_tests
mpirun -np 2 ./ex1
