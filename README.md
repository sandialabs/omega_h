![Omega_h Logo][1]

# Omega_h
> Reliable mesh adaptation

Omega_h is a C++11 library that implements tetrahedron and triangle mesh adaptativity,
with a focus on scalable HPC performance using (optionally) MPI, OpenMP, or CUDA.
It is intended to provided adaptive functionality to existing simulation codes.
Mesh adaptivity allows one to minimize both discretization error and number
of degrees of freedom live during the simulation, as well as enabling moving
object and evolving geometry simulations.
Omega_h will do this for you in a way that is fast, memory-efficient, and
portable across many different architectures.

## Installing / Getting started

For a bare minimum setup with no parallelism, you just need [CMake][0],
a C++11 compiler, and preferably [ZLib][6] installed.

```shell
git clone git@github.com:ibaned/omega_h2.git
cd omega_h2
cmake . -DCMAKE_INSTALL_PREFIX=/your/choice
make install
```

This should install Omega_h under the given prefix in a way you can
access from your own CMake files using these CMake commands:

```cmake
find_package(Omega_h)
target_link_libraries(myprogram Omega_h::omega_h)
```

## Features

Omega_h provides at least the following:
* Adaptation of tetrahedral and triangle meshes in parallel
* Anisotropic metric field support
* Given good input element quality, the output element
  quality is also guaranteed.
* Scalable MPI parallelism
* On-node OpenMP or CUDA parallelism using [Kokkos][2]
* Fully deterministic execution
* Given the same mesh, global numbering, and size field,
  results will be independent of parallel partitioning
  and ordering.

## Configuration

Below we document some key CMake configuration options:

#### Omega_h_USE_MPI
Default: `OFF`

Whether to enable MPI parallelism.
We recommend using [MPICH][3] or another MPI 3.0 implementation,
but we also support MPI version 2.1.
If this is `ON`, set `CMAKE_CXX_COMPILER` to your MPI compiler wrapper.

#### Omega_h_USE_Kokkos
Default: `OFF`

Whether to use [Kokkos][2] for on-node parallelism.
If this is `ON`, set `Kokkos_PREFIX` to your Kokkos installation.
You can install Kokkos as part of [Trilinos][4].
Please see [this file][8] for an example of how to configure
Trilinos to install only the Kokkos package.

#### Omega_h_USE_CUDA
Default: `OFF`

Whether to use CUDA via Kokkos.
This requires `Omega_h_USE_Kokkos`.
If this is `ON`, set `CMAKE_CXX_COMPILER` to your copy of
[nvcc_wrapper][7].

#### Omega_h_ONE_FILE
Default: `ON`

Omega_h includes all other sources in `omega_h.cpp` for a (up to 4X)
faster compile time from scratch.
Reasons to set this to `OFF` include:
* You have many cores on the compiling machine (at least 8), then
  you can get a faster compile time with `make -j`.
* You are a developer who would prefer faster incremental rebuilds.

## Contributing

Please open a Github issue to ask a question, report a bug,
request features, etc.
If you'd like to contribute, please fork the repository and use a feature
branch. Pull requests are welcome.

## Licensing

This library is released under the FreeBSD license.

[0]: https://cmake.org
[1]: https://raw.githubusercontent.com/ibaned/omega_h2/master/aux/omega_h.png
[2]: https://github.com/kokkos/kokkos
[3]: http://www.mpich.org
[4]: https://github.com/trilinos/Trilinos
[5]: http://clang.llvm.org/docs/ClangFormat.html
[6]: http://zlib.net
[7]: http://github.com/kokkos/nvcc_wrapper
[8]: https://github.com/ibaned/omega_h2/blob/master/aux/do-config-trilinos-kokkos.sh
