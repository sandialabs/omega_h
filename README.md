![omega_h logo][5] Omega_h
--------------------------

Omega_h is a C++11 library for HPC use of triangle and
tetrahedral meshes.
Its focus is on mesh adaptation: using a series of local
modifications, a mesh is altered to satisfy edge length and
element quality constraints while preserving the field
data stored on the mesh.

Omega_h can use the following forms of parallelism:

1. MPI distributed memory and processes
2. OpenMP shared memory and threads
3. CUDA shared memory and threads

The latter two are enabled via the [Kokkos library][1].
All forms of parallelism are optional, Omega_h can
run in serial without an MPI or Kokkos installation.

Omega_h has some very useful guarantees in terms of
reproducibility, which are important for regression
testing, debugging, and research publication:

1. Results are fully deterministic in serial and parallel,
even in OpenMP and CUDA modes.
2. Results are independent of the ordering of mesh entities
3. Results are independent of the partitioning of the mesh
4. Parallel behavior is independent of the local ordering
of mesh entities

Compiling
---------

The bare minimum requirements are the [CMake build system][2]
and a C++11 compliant compiler.
Our build system supports Clang, GNU, and Intel compilers.
It is strongly recommended to use the widely available [zlib][3]
for file compression.

If you use MPI, we recommend [MPICH][4] or another MPI 3.0
compliant implementation, although we do support the MPI 2.1
standard which still lingers on BlueGene/Q machines.

[1]: https://github.com/kokkos/kokkos
[2]: https://cmake.org
[3]: http://zlib.net
[4]: http://www.mpich.org
[5]: https://github.com/ibaned/omega_h2/raw/master/omega_h.png
