cmake_minimum_required(VERSION 3.0.0)
project(Omega_h CXX)

#trying to force BUILD_TESTING to default to OFF
set(BUILD_TESTING OFF CACHE STRING "Build and run tests")

include(CTest)
enable_testing()

set(OSH_MAJOR "1")
set(OSH_MINOR "5")
set(OSH_PATCH "0")

option(OSH_USE_MPI "Use MPI for parallelism" ON)
option(OSH_USE_KOKKOS "Use Kokkos for on-node parallelism" OFF)
option(OSH_USE_OPENMP "Use Kokkos+OpenMP for on-node parallelism" OFF)
option(OSH_USE_PTHREAD "Use Kokkos+Pthread for on-node parallelism" OFF)
option(OSH_USE_CUDA "Use Kokkos+CUDA for on-node parallelism" OFF)
option(OSH_USE_ZLIB "Use zlib for file compression" ON)
option(OSH_CHECK_BOUNDS "Check array bounds (makes code slow too)" OFF)
option(OSH_OPTIMIZE "Compile with optimization" ON)
option(OSH_SYMBOLS "Compile with debug symbols" ON)
option(OSH_BUILD_SMB2OSH "Enable building the smb2osh program" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(OSH_SANITIZE_ADDRESS "Use -fsanitize=address" OFF)
option(OSH_PROTECT "Catch OS signals and print stack" OFF)
set(Gmodel_PREFIX "" CACHE PATH "Gmodel installation directory")
set(Gmsh_PREFIX "" CACHE PATH "Gmsh installation directory")
option(OSH_ONE_FILE "Compile omega_h as a single file" ON)
set(OSH_CXX_FLAGS "" CACHE STRING "Extra compiler flags")
set(ZLIB_PREFIX "" CACHE PATH "zlib installation directory")

#If not building shared libs, then prefer static
#dependency libs
if(NOT BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".dylib")
endif()

set(FLAGS "--std=c++11")
if(OSH_OPTIMIZE)
  set(FLAGS "${FLAGS} -O2")
else()
  set(FLAGS "${FLAGS} -O0")
endif()
if(OSH_SYMBOLS)
  set(FLAGS "${FLAGS} -g")
endif()
if(OSH_SANITIZE_ADDRESS)
  set(FLAGS "${FLAGS} -fsanitize=address")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
  set(FLAGS "${FLAGS} -Werror -Weverything -Wno-c++98-compat")
  set(FLAGS "${FLAGS} -Wno-c++98-compat-pedantic -Wno-padded")
  set(FLAGS "${FLAGS} -Wno-float-equal")
  set(FLAGS "${FLAGS} -Wno-weak-template-vtables")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if(OSH_USE_CUDA)
    set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
    set(FLAGS "${FLAGS} -expt-extended-lambda")
  elseif(OSH_USE_OPENMP)
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
    set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
    set(FLAGS "${FLAGS} -fopenmp")
  elseif(OSH_USE_PTHREAD)
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
    set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
    set(FLAGS "${FLAGS} -pthread")
  else()
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer")
    set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
  endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(FLAGS "${FLAGS} -debug")
  if(OSH_USE_OPENMP)
    set(FLAGS "${FLAGS} -qopenmp")
  endif()
else()
  message(FATAL_ERROR "Unexpected compiler type ${CMAKE_CXX_COMPILER_ID}")
endif()
set(CMAKE_CXX_FLAGS "${FLAGS} ${OSH_CXX_FLAGS}")
message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")

# CMake RPATH "always full" configuration, see:
# https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH False)
# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH False)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)

add_subdirectory(src)

configure_file(
  "${PROJECT_SOURCE_DIR}/omega_h-config.cmake.in"
  "${PROJECT_BINARY_DIR}/omega_h-config.cmake"
  @ONLY)
configure_file(
  "${PROJECT_SOURCE_DIR}/omega_h-config-version.cmake.in"
  "${PROJECT_BINARY_DIR}/omega_h-config-version.cmake"
  @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}/omega_h-config.cmake"
  "${PROJECT_BINARY_DIR}/omega_h-config-version.cmake"
  DESTINATION lib/omega_h)
