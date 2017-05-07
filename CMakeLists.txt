cmake_minimum_required(VERSION 3.0.0)

project(Omega_h VERSION 9.0.0 LANGUAGES NONE)

set(USE_XSDK_DEFAULTS_DEFAULT OFF)

include(cmake/XSDKDefaults.cmake)

enable_language(CXX)

include(cmake/bob.cmake)

bob_begin_package()

option(Omega_h_USE_MPI "Use MPI for parallelism" OFF)
message(STATUS "Omega_h_USE_MPI: ${Omega_h_USE_MPI}")
option(Omega_h_USE_OpenMP "Use Kokkos+OpenMP for on-node parallelism" OFF)
option(Omega_h_USE_PTHREADS "Use Kokkos+Pthread for on-node parallelism" OFF)
option(Omega_h_USE_CUDA "Use Kokkos+CUDA for on-node parallelism" OFF)
option(Omega_h_CHECK_BOUNDS "Check array bounds (makes code slow too)" OFF)
message(STATUS "Omega_h_CHECK_BOUNDS: ${Omega_h_CHECK_BOUNDS}")
option(Omega_h_SANITIZE_ADDRESS "Use -fsanitize=address" OFF)
option(Omega_h_PROTECT "Catch OS signals and print stack" OFF)
set(Gmsh_PREFIX "" CACHE PATH "Gmsh installation directory")
option(Omega_h_ONE_FILE "Compile Omega_h as a single file" OFF)
message(STATUS "Omega_h_ONE_FILE: ${Omega_h_ONE_FILE}")
option(Omega_h_USE_EGADS "Use EGADS from ESP for geometry" OFF)
set(EGADS_PREFIX "" CACHE PATH "EGADS (or ESP) installation directory")
set(Omega_h_DATA "" CACHE PATH "Path to omega_h-data test files")
set(Trilinos_PREFIX "" CACHE PATH "Path to Trilinos install")

bob_begin_cxx_flags()
bob_cxx11_flags()
if(Omega_h_SANITIZE_ADDRESS)
  set(FLAGS "${FLAGS} -fsanitize=address")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if(Omega_h_USE_CUDA)
    set(FLAGS "${FLAGS} -expt-extended-lambda")
  elseif(Omega_h_USE_OpenMP)
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
    set(FLAGS "${FLAGS} -fopenmp")
  elseif(Omega_h_USE_PTHREADS)
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
    set(FLAGS "${FLAGS} -pthread")
  else()
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
  endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
  if(Omega_h_USE_OpenMP)
    set(FLAGS "${FLAGS} -qopenmp")
  endif()
else()
  message(WARNING "Unexpected compiler type ${CMAKE_CXX_COMPILER_ID}")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
bob_end_cxx_flags()

set(Omega_h_USE_ZLIB_DEFAULT ON)
bob_public_dep(ZLIB)

set(Omega_h_USE_Kokkos_DEFAULT OFF)
set(Kokkos_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(Kokkos)

if(Omega_h_USE_Kokkos)
  get_target_property(kokkoscore_LOCATION kokkoscore LOCATION)
  message(STATUS "kokkoscore LOCATION: ${kokkoscore_LOCATION}")
  get_target_property(kokkoscore_LINK kokkoscore INTERFACE_LINK_LIBRARIES)
  if(kokkoscore_LINK)
    message(STATUS "kokkoscore INTERFACE_LINK_LIBRARIES: ${kokkoscore_LINK}")
  endif()
endif()

set(libMeshb_REQUIRED_VERSION 7.24)
set(Omega_h_USE_libMeshb_DEFAULT OFF)
bob_public_dep(libMeshb)

set(Omega_h_USE_SEACASExodus_DEFAULT OFF)
set(SEACASExodus_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(SEACASExodus)

set(Omega_h_USE_TeuchosParameterList_DEFAULT OFF)
set(TeuchosParameterList_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(TeuchosParameterList)

set(Omega_h_USE_TeuchosComm_DEFAULT OFF)
set(TeuchosComm_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(TeuchosComm)

if(Omega_h_USE_TeuchosParameterList AND Omega_h_USE_TeuchosComm)
  message(STATUS "Using Teuchos")
  set(USE_TEUCHOS ON)
elseif(Omega_h_USE_TeuchosParameterList OR Omega_h_USE_TeuchosComm)
  message(FATAL_ERROR "Must use both TeuchosParameterList and TeuchosComm or neither")
else()
  message(STATUS "Not using Teuchos")
  set(USE_TEUCHOS OFF)
endif()

set(Gmodel_REQUIRED_VERSION 2.1.0)
set(Omega_h_USE_Gmodel_DEFAULT OFF)
bob_private_dep(Gmodel)

set(Omega_h_KEY_BOOLS
    Omega_h_USE_MPI
    Omega_h_USE_Kokkos
    Omega_h_USE_OpenMP
    Omega_h_USE_CUDA
    Omega_h_USE_ZLIB
    Omega_h_USE_libMeshb
    Omega_h_USE_EGADS
    Omega_h_USE_SEACASExodus
    Omega_h_CHECK_BOUNDS
    Omega_h_PROTECT
   )

set(Omega_h_KEY_INTS
    Omega_h_VERSION_MAJOR
    Omega_h_VERSION_MINOR
    Omega_h_VERSION_PATCH
   )

set(Omega_h_KEY_STRINGS
    Omega_h_SEMVER
   )

bob_form_semver()

add_subdirectory(src)

bob_end_package()
