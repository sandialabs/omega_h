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
git clone git@github.com:ibaned/omega_h.git
cd omega_h
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

Omega\_h provides at least the following:
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

#### Omega_h_USE_Trilinos
Default: `OFF`

Whether to use the Kokkos and Teuchos Trilinos packages.
If this is `ON`, set `Trilinos_PREFIX` to your Trilinos installation.
Kokkos provides on-node parallelism, and Omeg\_h will use the same
default execution space that Kokkos was built with.
Teuchos provides parameter lists and file I/O for them,
which are usable through `Omega_h_teuchos.hpp`.

#### Omega_h_USE_SEACASExodus
Default: `OFF`

Whether to use the Exodus subpackage of the SEACAS Trilinos package.
This allows reading and writing Exodus files from Omega\_h.
By default, it will look for this dependency in `Trilinos_PREFIX`.

## Contributing

Please open a Github issue to ask a question, report a bug,
request features, etc.
If you'd like to contribute, please fork the repository and use a feature
branch. Pull requests are welcome.

## Licensing

This library is released under the FreeBSD license.

[0]: https://cmake.org
[1]: https://raw.githubusercontent.com/ibaned/omega_h/master/aux/omega_h.png
[2]: https://github.com/kokkos/kokkos
[3]: http://www.mpich.org
[4]: https://github.com/trilinos/Trilinos
[5]: http://clang.llvm.org/docs/ClangFormat.html
[6]: http://zlib.net
[7]: http://github.com/kokkos/nvcc_wrapper
[8]: https://github.com/ibaned/omega_h/blob/master/aux/do-config-trilinos-kokkos.sh
