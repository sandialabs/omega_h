cmake_minimum_required(VERSION 3.18.0...${CMAKE_VERSION})

project(Omega_h VERSION 10.5.0 LANGUAGES CXX)
set(Omega_h_USE_STK_DEFAULT OFF)
option(Omega_h_USE_STK "Whether to build the STK interface" ${Omega_h_USE_STK_DEFAULT})

set(Omega_h_USE_DOLFIN_DEFAULT OFF)
option(Omega_h_USE_DOLFIN "Whether to build the DOLFIN interface" ${Omega_h_USE_DOLFIN_DEFAULT})
if (Omega_h_USE_DOLFIN)
  # The FindPETSc module will use try_compile with a C code file
  enable_language(C)
endif()

include(cmake/bob.cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

bob_begin_package()

bob_option(Omega_h_CHECK_BOUNDS "Check array bounds when running on host (makes code slow too)" OFF)
bob_option(Omega_h_THROW "Errors throw exceptions instead of abort" ${USE_XSDK_DEFAULTS})
bob_input(Omega_h_DATA "" PATH "Path to omega_h-data test files")
bob_option(Omega_h_USE_EGADS "Use EGADS from ESP for geometry" OFF)
bob_input(EGADS_PREFIX "" PATH "EGADS (or ESP) installation directory")
bob_option(Omega_h_USE_Kokkos "Use Kokkos as a backend" OFF)
bob_input(Kokkos_PREFIX "" PATH "Path to Kokkos install")
bob_option(Omega_h_USE_CUDA_AWARE_MPI "Assume MPI is CUDA-aware, make use of that" OFF)
bob_option(Omega_h_USE_SimModSuite "Enable reading Simmetrix MeshSim meshes" OFF)
bob_input(Omega_h_VALGRIND "" STRING "Valgrind plus arguments for testing")
bob_option(Omega_h_EXAMPLES "Compile examples" OFF)

bob_option(Omega_h_USE_MPI "Use MPI for parallelism" OFF)
bob_option(Omega_h_ENABLE_DEMANGLED_STACKTRACE "Add linker options to enable human-readable stacktraces on gnu-c platforms." OFF)
bob_option(Omega_h_DBG "Enable debug prints, stacktraces, etc." OFF)
bob_option(Omega_h_GPU_CHECK "Run GPU check after each test" OFF) #the check command is hardcoded!

if (Omega_h_ENABLE_DEMANGLED_STACKTRACE)
  message(STATUS "CMAKE_BUILD_TYPE= ${CMAKE_BUILD_TYPE}")
  if ((CMAKE_BUILD_TYPE STREQUAL "Debug") OR (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo"))
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic -ldl")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
    message(STATUS "Omega_h adding -rdynamic to CMAKE_EXE_LINKER_FLAGS= ${CMAKE_EXE_LINKER_FLAGS}")
  endif()
endif()

if (Omega_h_USE_MPI)
  enable_language(C)
  set(MPI_REQUIRED_VERSION 3)
  bob_add_dependency(PUBLIC NAME MPI COMPONENTS CXX TARGETS MPI::MPI_CXX)
endif()

set(Omega_h_USE_ZLIB_DEFAULT ON)
bob_add_dependency(PUBLIC NAME ZLIB TARGETS ZLIB::ZLIB)

set(Omega_h_USE_Kokkos_DEFAULT OFF)
set(KokkosCore_PREFIX_DEFAULT ${Kokkos_PREFIX})
bob_add_dependency(PUBLIC NAME Kokkos TARGETS Kokkos::kokkos)
if (Omega_h_USE_Kokkos)
  set_target_properties(Kokkos::kokkos PROPERTIES INTERFACE_COMPILE_OPTIONS "")
endif()

set(libMeshb_REQUIRED_VERSION 7.24)
set(Omega_h_USE_libMeshb_DEFAULT OFF)
bob_add_dependency(PUBLIC NAME libMeshb TARGETS libMeshb::Meshb.7)

set(Omega_h_USE_Gmsh_DEFAULT OFF)
bob_add_dependency(PUBLIC NAME Gmsh
    INCLUDE_DIR_VARS Gmsh_INCLUDE_DIRS
    LIBRARY_VARS Gmsh_LIBRARIES)
if(Omega_h_USE_Gmsh)
  if(Gmsh_VERSION_STRING VERSION_LESS 4.4.1)
    message(SEND_ERROR "Gmsh minimum version is 4.4.1 but found version ${Gmsh_VERSION_STRING}")
  endif()
endif()

set(Gmodel_REQUIRED_VERSION 2.1.0)
set(Omega_h_USE_Gmodel_DEFAULT OFF)
bob_add_dependency(PRIVATE NAME Gmodel TARGETS gmodel)

set(Omega_h_USE_SEACASExodus_DEFAULT OFF)
set(SEACASExodus_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_add_dependency(PUBLIC NAME SEACASExodus TARGETS SEACASExodus::all_libs)

set(Omega_h_USE_pybind11_DEFAULT OFF)
bob_public_dep(pybind11)

if(Omega_h_USE_Kokkos)
  if ("CUDA" IN_LIST Kokkos_DEVICES)
      add_definitions(-DOmega_h_USE_CUDA)
      bob_option(Omega_h_USE_CUDA "Whether to use CUDA" "ON")
  endif()
  if ("OPENMP" IN_LIST Kokkos_DEVICES)
      add_definitions(-DOmega_h_USE_OpenMP)
      bob_option(Omega_h_USE_OpenMP "Whether to use OpenMP" "ON")
  endif()
  if ("SYCL" IN_LIST Kokkos_DEVICES)
      add_definitions(-DOmega_h_USE_SYCL)
      bob_option(Omega_h_USE_SYCL "Whether to use SYCL" "ON")
  endif()
  if ("HIP" IN_LIST Kokkos_DEVICES)
    add_definitions(-DOmega_h_USE_HIP)
    bob_option(Omega_h_USE_HIP "Whether to use HIP" "ON")
    endif()
  if (Omega_h_USE_CUDA AND (NOT "CUDA_LAMBDA" IN_LIST Kokkos_OPTIONS))
    message(FATAL_ERROR
            "Please reconfigure Kokkos with -DKokkos_ENABLE_Cuda_Lambda:BOOL=ON")
  endif()
endif()

if (Omega_h_USE_KOKKOS_MEM_SPACE STREQUAL "DEVICE")
  add_definitions(-DOmega_h_MEM_SPACE_DEVICE)
endif()

bob_option(Omega_h_USE_MPI "Use MPI for parallelism" OFF)
if(Omega_h_USE_SimModSuite AND NOT Omega_h_USE_MPI)
  message(FATAL_ERROR "Omega_h_USE_MPI must be enabled if Omega_h_USE_SimModSuite is enabled")
endif()
bob_option(Omega_h_USE_OpenMP "Whether to use OpenMP" "${Kokkos_HAS_OpenMP}")
bob_option(Omega_h_USE_CUDA "Whether to use CUDA" "${Kokkos_HAS_CUDA}")

if (Omega_h_USE_CUDA)
  enable_language(CUDA)
  #check for versions of cuda/thrust that dont't work with omegah
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.2
      AND CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 11.3)
    message(FATAL_ERROR "CUDA 11.2 does not support Omega_h, "
      "use an older or newer version")
  endif()
  #check for bad combos
  if(Omega_h_CHECK_BOUNDS)
    message(FATAL_ERROR "Bounds checking is only supported on the host, "
      "disable CUDA/Kokkos+CUDA or bounds checking")
  endif()
endif()

include(cmake/osh_use_dolfin.cmake)
osh_use_dolfin()

set(Omega_h_IS_SHARED ${BUILD_SHARED_LIBS})

set(Omega_h_KEY_BOOLS
    Omega_h_USE_MPI
    Omega_h_ENABLE_DEMANGLED_STACKTRACE
    Omega_h_DBG
    Omega_h_USE_Kokkos
    Omega_h_USE_OpenMP
    Omega_h_USE_CUDA
    Omega_h_USE_SYCL
    Omega_h_USE_HIP
    Omega_h_MEM_SPACE_DEVICE
    Omega_h_USE_ZLIB
    Omega_h_USE_libMeshb
    Omega_h_USE_EGADS
    Omega_h_USE_SEACASExodus
    Omega_h_USE_SimModSuite
    Omega_h_USE_DOLFIN
    Omega_h_USE_dwarf
    Omega_h_CHECK_BOUNDS
    Omega_h_THROW
    Omega_h_USE_CUDA_AWARE_MPI
    Omega_h_USE_Gmsh
    Omega_h_IS_SHARED
   )

set(Omega_h_KEY_INTS
    Omega_h_VERSION_MAJOR
    Omega_h_VERSION_MINOR
    Omega_h_VERSION_PATCH
   )

set(Omega_h_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(Omega_h_KEY_STRINGS
    Omega_h_SEMVER
    Omega_h_COMMIT
    Omega_h_CXX_FLAGS
    Omega_h_CMAKE_ARGS
   )

bob_get_commit()
bob_get_semver()

if(Omega_h_EXAMPLES)
  find_program(GMSH_EXECUTABLE gmsh REQUIRED)
  add_subdirectory(example/field_on_square)
  add_subdirectory(example/castle)
  add_subdirectory(example/gradient)
  add_subdirectory(example/laplacian)
  if(Omega_h_USE_DOLFIN)
    set(SAVED_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    add_subdirectory(example/fenics/cpp)
    set(CMAKE_CXX_FLAGS "${SAVED_CXX_FLAGS}")
  endif()
endif()
add_subdirectory(src)

if (Omega_h_USE_CUDA)
  install(PROGRAMS
      ${CMAKE_CURRENT_SOURCE_DIR}/bin/nvcc_wrapper
      DESTINATION bin)
endif()


bob_end_package()
