cmake_minimum_required(VERSION 3.0.0)

project(Omega_h VERSION 9.13.13 LANGUAGES CXX)

set(Omega_h_USE_DOLFIN_DEFAULT OFF)
option(Omega_h_USE_DOLFIN "Whether to build the DOLFIN interface" ${Omega_h_USE_DOLFIN_DEFAULT})
if (Omega_h_USE_DOLFIN)
  # The FindPETSc module will use try_compile with a C code file
  enable_language(C)
endif()

include(cmake/bob.cmake)

bob_begin_package()

bob_option(Omega_h_CHECK_BOUNDS "Check array bounds (makes code slow too)" OFF)
bob_option(Omega_h_THROW "Errors throw exceptions instead of abort" ${USE_XSDK_DEFAULTS})
bob_input(Omega_h_DATA "" PATH "Path to omega_h-data test files")
bob_option(Omega_h_USE_EGADS "Use EGADS from ESP for geometry" OFF)
bob_input(EGADS_PREFIX "" PATH "EGADS (or ESP) installation directory")
bob_option(Omega_h_USE_Trilinos "Use Teuchos and Kokkos" OFF)
bob_input(Trilinos_PREFIX "" PATH "Path to Trilinos install")
bob_option(Omega_h_USE_CUDA_AWARE_MPI "Assume MPI is CUDA-aware, make use of that" OFF)
bob_input(Omega_h_VALGRIND "" STRING "Valgrind plus arguments for testing")
bob_option(Omega_h_EXAMPLES "Compile examples" OFF)

set(Omega_h_USE_ZLIB_DEFAULT ON)
bob_add_dependency(PUBLIC NAME ZLIB TARGETS ZLIB::ZLIB)

set(Omega_h_USE_KokkosCore_DEFAULT ${Omega_h_USE_Trilinos})
set(KokkosCore_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_add_dependency(PUBLIC NAME KokkosCore TARGETS kokkoscore INCLUDE_DIR_VARS KokkosCore_INCLUDE_DIRS)

set(libMeshb_REQUIRED_VERSION 7.24)
set(Omega_h_USE_libMeshb_DEFAULT OFF)
bob_add_dependency(PUBLIC NAME libMeshb TARGETS libMeshb::Meshb.7)

set(Gmodel_REQUIRED_VERSION 2.1.0)
set(Omega_h_USE_Gmodel_DEFAULT OFF)
bob_add_dependency(PRIVATE NAME Gmodel TARGETS gmodel)

set(Omega_h_USE_SEACASExodus_DEFAULT OFF)
set(SEACASExodus_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_add_dependency(PUBLIC NAME SEACASExodus
    TARGETS exodus
    INCLUDE_DIR_VARS
      SEACASExodus_INCLUDE_DIRS
      SEACASExodus_TPL_INCLUDE_DIRS
    LIBRARY_VARS SEACASExodus_TPL_LIBRARIES)

set(Omega_h_USE_TeuchosParameterList_DEFAULT ${Omega_h_USE_Trilinos})
set(TeuchosParameterList_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_add_dependency(PUBLIC NAME TeuchosParameterList
    TARGETS teuchosparameterlist
    INCLUDE_DIR_VARS
      TeuchosParameterList_INCLUDE_DIRS
      TeuchosParameterList_TPL_INCLUDE_DIRS
    LIBRARY_VARS TeuchosParameterList_TPL_LIBRARIES)

set(Omega_h_USE_TeuchosComm_DEFAULT ${Omega_h_USE_Trilinos})
set(TeuchosComm_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_add_dependency(PUBLIC NAME TeuchosComm
    TARGETS teuchoscomm
    INCLUDE_DIR_VARS
      TeuchosComm_INCLUDE_DIRS
      TeuchosComm_TPL_INCLUDE_DIRS
    LIBRARY_VARS TeuchosComm_TPL_LIBRARIES)

set(Omega_h_USE_TeuchosParser_DEFAULT ${Omega_h_USE_Trilinos})
set(TeuchosParser_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_add_dependency(PUBLIC NAME TeuchosParser TARGETS teuchosparser
    INCLUDE_DIR_VARS
      TeuchosParser_INCLUDE_DIRS
    )

set(Omega_h_USE_pybind11_DEFAULT OFF)
bob_public_dep(pybind11)

bob_option(Omega_h_USE_dwarf "Use libdwarf for better stack traces" OFF)
if (Omega_h_USE_dwarf)
	find_path(LIBDWARF_INCLUDE_DIR NAMES "libdwarf.h" PATH_SUFFIXES libdwarf)
	find_path(LIBELF_INCLUDE_DIR NAMES "libelf.h")
	find_path(LIBDL_INCLUDE_DIR NAMES "dlfcn.h")
	find_library(LIBDWARF_LIBRARY dwarf)
	find_library(LIBELF_LIBRARY elf)
	find_library(LIBDL_LIBRARY dl)
	set(LIBDWARF_INCLUDE_DIRS ${LIBDWARF_INCLUDE_DIR} ${LIBELF_INCLUDE_DIR} ${LIBDL_INCLUDE_DIR})
	set(LIBDWARF_LIBRARIES ${LIBDWARF_LIBRARY} ${LIBELF_LIBRARY} ${LIBDL_LIBRARY})
endif()

if  ( Omega_h_USE_TeuchosParameterList
  AND Omega_h_USE_TeuchosComm
  AND Omega_h_USE_TeuchosParser
  )
  message(STATUS "Using Teuchos")
  set(Omega_h_USE_Teuchos ON)
elseif  ( Omega_h_USE_TeuchosParameterList
       OR Omega_h_USE_TeuchosComm
       OR Omega_h_USE_TeuchosParser
       )
  message(FATAL_ERROR "Must use TeuchosParameterList, TeuchosComm, and TeuchosParser or none")
else()
  message(STATUS "Not using Teuchos")
  set(Omega_h_USE_Teuchos OFF)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/cmake/detect_trilinos_opts.cmake)
detect_trilinos_opts()
if (KokkosCore_HAS_CUDA AND (NOT KokkosCore_HAS_CUDA_LAMBDA))
  message(FATAL_ERROR
          "Please reconfigure Kokkos with -DKokkos_ENABLE_Cuda_Lambda:BOOL=ON")
endif()

bob_option(Omega_h_USE_MPI "Use MPI for parallelism" "${TeuchosComm_HAS_MPI}")
bob_option(Omega_h_USE_OpenMP "Whether to use OpenMP" "${KokkosCore_HAS_OpenMP}")
bob_option(Omega_h_USE_CUDA "Whether to use CUDA" "${KokkosCore_HAS_CUDA}")

if (Omega_h_USE_CUDA)
  find_package(CUDA)
  bob_input(Omega_h_CUDA_ARCH "" STRING "CUDA arch specifier")
endif()

bob_begin_cxx_flags()
bob_cxx11_flags()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
  if(Omega_h_USE_CUDA)
    if (Omega_h_CUDA_ARCH)
      add_compile_options($<$<COMPILE_LANGUAGE:CXX>:--cuda-gpu-arch=${Omega_h_CUDA_ARCH}>)
    endif()
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-x>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:cuda>)
    # -Wweak-vtables gives false positives with -x cuda
    set(FLAGS "${FLAGS} -Wno-weak-vtables")
  endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if(Omega_h_USE_CUDA)
    set(FLAGS "${FLAGS} -expt-extended-lambda -lineinfo")
    if (Omega_h_CUDA_ARCH)
      set(FLAGS "${FLAGS} -arch=${Omega_h_CUDA_ARCH}")
    endif()
    if (CUDA_VERSION EQUAL 9)
      set(FLAGS "${FLAGS} -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored")
    endif()
  else()
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
    if(Omega_h_USE_OpenMP)
      set(FLAGS "${FLAGS} -fopenmp")
    endif()
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

include(cmake/osh_use_dolfin.cmake)
osh_use_dolfin()

set(Omega_h_KEY_BOOLS
    Omega_h_USE_MPI
    Omega_h_USE_KokkosCore
    Omega_h_USE_OpenMP
    Omega_h_USE_CUDA
    Omega_h_USE_ZLIB
    Omega_h_USE_libMeshb
    Omega_h_USE_EGADS
    Omega_h_USE_SEACASExodus
    Omega_h_USE_Teuchos
    Omega_h_USE_DOLFIN
    Omega_h_USE_dwarf
    Omega_h_CHECK_BOUNDS
    Omega_h_THROW
    Omega_h_USE_CUDA_AWARE_MPI
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
endif()
add_subdirectory(src)

bob_end_package()
