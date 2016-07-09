![Omega_h Logo][1]

# Omega_h
> Reliable mesh adaptation

Omega_h is a C++11 library that implements tet and triangle mesh adaptativity,
with a focus on scalable HPC performance using MPI and OpenMP or CUDA.
It is intended to provided adaptive functionality to existing simulation codes.
Mesh adaptivity allows one to minimize both discretization error and number
of degrees of freedom live during the simulation, as well as enabling moving
object and evolving geometry simulations.
Omega_h will do this for you in a way that is fast, memory-efficient, and
portable across many different architectures.

## Installing / Getting started

For a bare minimum setup with no parallelism, you just need [CMake][0],
a C++11 compiler, and [ZLib][6] installed.

```shell
git clone git@github.com:ibaned/omega_h2.git
cd omega_h2
cmake . -DCMAKE_INSTALL_PREFIX=/your/choice -DOSH_USE_MPI=OFF
make install
```

This should install Omega_h under the given prefix in a way you can
access from your own CMake files using these CMake commands:

```cmake
find_package(Omega_h 1.0.0)
include_directories(${Omega_h_INCLUDE_DIRS})
target_link_libraries(myprogram omega_h)
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

#### OSH_USE_MPI
Default: `ON`

Whether to enable MPI parallelism.
We recommend using [MPICH][3] or another MPI 3.0 implementation,
but we also support MPI version 2.1.
If this is `ON`, set `CMAKE_CXX_COMPILER` to your MPI compiler wrapper.

#### OSH_USE_KOKKOS
Default: `OFF`

Whether to use [Kokkos][2] for on-node parallelism.
If this is `ON`, set `Kokkos_PREFIX` to your Kokkos installation.
You can install Kokkos either on [its own][2] or as part of
[Trilinos][4].

#### OSH_USE_CUDA
Default: `OFF`

Whether to use CUDA via Kokkos.
This requires `OSH_USE_KOKKOS`.
If this is `ON`, set `OSH_CUDA_ARCH` to the compute capability of
your GPUs.

## Contributing

Please open a Github issue to ask a question, report a bug,
request features, etc.
If you'd like to contribute, please fork the repository and use a feature
branch. Pull requests are welcome.

We use [clang-format][5] to maintain a consistent code style.
Developers are encouraged to install it

```shell
sudo apt-get install clang-format-3.5
```

And run something like this before committing unstaged changes:

```shell
git clang-format -f
```

## Licensing

This library is released under the FreeBSD license.

[0]: https://cmake.org
[1]: https://raw.githubusercontent.com/ibaned/omega_h2/master/omega_h.png
[2]: https://github.com/kokkos/kokkos
[3]: http://www.mpich.org
[4]: https://github.com/trilinos/Trilinos
[5]: http://clang.llvm.org/docs/ClangFormat.html
[6]: http://zlib.net
