cmake_minimum_required(VERSION 3.0.0)

project(Omega_h VERSION 9.1.2 LANGUAGES NONE)

set(USE_XSDK_DEFAULTS_DEFAULT OFF)

include(cmake/XSDKDefaults.cmake)

enable_language(CXX)

include(cmake/bob.cmake)

bob_begin_package()

option(Omega_h_CHECK_BOUNDS "Check array bounds (makes code slow too)" OFF)
message(STATUS "Omega_h_CHECK_BOUNDS: ${Omega_h_CHECK_BOUNDS}")
option(Omega_h_ONE_FILE "Compile Omega_h as a single file" OFF)
message(STATUS "Omega_h_ONE_FILE: ${Omega_h_ONE_FILE}")
set(Omega_h_DATA "" CACHE PATH "Path to omega_h-data test files")
option(Omega_h_USE_EGADS "Use EGADS from ESP for geometry" OFF)
set(EGADS_PREFIX "" CACHE PATH "EGADS (or ESP) installation directory")
option(Omega_h_USE_Trilinos "Use Teuchos and Kokkos" OFF)
message(STATUS "Omega_h_USE_Trilinos: ${Omega_h_USE_Trilinos}")
set(Trilinos_PREFIX "" CACHE PATH "Path to Trilinos install")

set(Omega_h_USE_ZLIB_DEFAULT ON)
bob_public_dep(ZLIB)

set(Omega_h_USE_KokkosCore_DEFAULT ${Omega_h_USE_Trilinos})
set(KokkosCore_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(KokkosCore)

if(Omega_h_USE_KokkosCore)
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

set(Gmodel_REQUIRED_VERSION 2.1.0)
set(Omega_h_USE_Gmodel_DEFAULT OFF)
bob_private_dep(Gmodel)

set(Omega_h_USE_SEACASExodus_DEFAULT OFF)
set(SEACASExodus_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(SEACASExodus)

set(Omega_h_USE_TeuchosParameterList_DEFAULT ${Omega_h_USE_Trilinos})
set(TeuchosParameterList_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(TeuchosParameterList)

set(Omega_h_USE_TeuchosComm_DEFAULT ${Omega_h_USE_Trilinos})
set(TeuchosComm_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(TeuchosComm)

set(Omega_h_USE_TeuchosParser_DEFAULT ${Omega_h_USE_Trilinos})
set(TeuchosParser_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_public_dep(TeuchosParser)

if  ( Omega_h_USE_TeuchosParameterList
  AND Omega_h_USE_TeuchosComm
  AND Omega_h_USE_TeuchosParser
  )
  message(STATUS "Using Teuchos")
  set(USE_TEUCHOS ON)
elseif  ( Omega_h_USE_TeuchosParameterList 
       OR Omega_h_USE_TeuchosComm
       OR Omega_h_USE_TeuchosParser
       )
  message(FATAL_ERROR "Must use TeuchosParameterList, TeuchosComm, and TeuchosParser or none")
else()
  message(STATUS "Not using Teuchos")
  set(USE_TEUCHOS OFF)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmake/detect_trilinos_opts.cmake)
detect_trilinos_opts()
if (KokkosCore_HAS_CUDA AND (NOT KokkosCore_HAS_CUDA_LAMBDA))
  message(FATAL_ERROR
          "Please reconfigure Kokkos with -DKokkos_ENABLE_Cuda_Lambda:BOOL=ON")
endif()

set(Omega_h_USE_MPI_DEFAULT OFF)
if (TeuchosComm_HAS_MPI)
  set(Omega_h_USE_MPI_DEFAULT ON)
endif()

option(Omega_h_USE_MPI "Use MPI for parallelism" ${Omega_h_USE_MPI_DEFAULT})
message(STATUS "Omega_h_USE_MPI: ${Omega_h_USE_MPI}")
set(Omega_h_USE_OpenMP ${KokkosCore_HAS_OpenMP})
message(STATUS "Omega_h_USE_OpenMP: ${Omega_h_USE_OpenMP}")
set(Omega_h_USE_CUDA ${KokkosCore_HAS_CUDA})
message(STATUS "Omega_h_USE_CUDA: ${Omega_h_USE_CUDA}")
set(Omega_h_USE_YAML ${TeuchosParameterList_HAS_YAML})
message(STATUS "Omega_h_USE_YAML: ${Omega_h_USE_YAML}")

bob_begin_cxx_flags()
bob_cxx11_flags()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if(Omega_h_USE_CUDA)
    set(FLAGS "${FLAGS} -expt-extended-lambda")
  elseif(Omega_h_USE_OpenMP)
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
    set(FLAGS "${FLAGS} -fopenmp")
  else()
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
  endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
  if(Omega_h_USE_OpenMP)
    set(FLAGS "${FLAGS} -fopenmp")
  endif()
else()
  message(WARNING "Unexpected compiler type ${CMAKE_CXX_COMPILER_ID}")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
bob_end_cxx_flags()

set(Omega_h_KEY_BOOLS
    Omega_h_USE_MPI
    Omega_h_USE_KokkosCore
    Omega_h_USE_OpenMP
    Omega_h_USE_CUDA
    Omega_h_USE_ZLIB
    Omega_h_USE_libMeshb
    Omega_h_USE_EGADS
    Omega_h_USE_SEACASExodus
    Omega_h_USE_TeuchosParameterList
    Omega_h_USE_TeuchosComm
    Omega_h_USE_TeuchosParser
    Omega_h_USE_YAML
    Omega_h_CHECK_BOUNDS
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
