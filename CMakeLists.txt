cmake_minimum_required(VERSION 2.6)
project(Omega_h CXX)

set(OSH_MAJOR "0")
set(OSH_MINOR "5")
set(OSH_PATCH "0")

option(OSH_USE_MPI "Use MPI for parallelism" ON)
option(OSH_USE_KOKKOS "Use Kokkos for on-node parallelism" OFF)
option(OSH_USE_OPENMP "Use Kokkos+OpenMP for on-node parallelism" OFF)
option(OSH_USE_CUDA "Use Kokkos+CUDA for on-node parallelism" OFF)
set(OSH_CUDA_ARCH "37" CACHE STRING "compute capability for NVCC -arch=sm_?")
option(OSH_USE_ZLIB "Use zlib for file compression" ON)
option(OSH_CHECK_BOUNDS "Check array bounds (makes code slow too)" OFF)
option(OSH_OPTIMIZE "Compile with optimization" ON)
option(OSH_SYMBOLS "Compile with debug symbols" ON)

if(OSH_USE_KOKKOS)
  find_path(KOKKOS_INCLUDE Kokkos_Core.hpp)
  find_library(KOKKOS_LIB NAMES kokkos kokkoscore)
  set(DEP_INCLUDE_DIRS ${DEP_INCLUDE_DIRS} ${KOKKOS_INCLUDE})
endif()

if(OSH_USE_ZLIB)
  find_package(ZLIB REQUIRED)
  set(DEP_INCLUDE_DIRS ${DEP_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
endif()

include_directories(${DEP_INCLUDE_DIRS})

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
set(FLAGS "--std=c++11 -fno-omit-frame-pointer -ferror-limit=1")
set(FLAGS "${FLAGS} -Werror -Weverything -Wno-c++98-compat")
set(FLAGS "${FLAGS} -Wno-c++98-compat-pedantic -Wno-padded")
set(FLAGS "${FLAGS} -Wno-float-equal")
set(FLAGS "${FLAGS} -Wno-weak-template-vtables")
if(OSH_OPTIMIZE)
  set(FLAGS "${FLAGS} -O2")
else()
  set(FLAGS "${FLAGS} -O0")
endif()
if(OSH_SYMBOLS)
  set(FLAGS "${FLAGS} -g")
endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
if(OSH_USE_CUDA)
set(FLAGS "--std=c++11 -O2 -g -Werror -Wall -Wextra")
set(FLAGS "${FLAGS} -expt-extended-lambda -arch=sm_${OSH_CUDA_ARCH}")
elseif(OSH_USE_OPENMP)
set(FLAGS "--std=c++11 -O2 -g -fno-omit-frame-pointer -fmax-errors=1")
set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
set(FLAGS "${FLAGS} -fopenmp")
else()
set(FLAGS "--std=c++11 -O2 -g -fno-omit-frame-pointer -fmax-errors=1")
set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
endif()
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
set(FLAGS "-std=c++11 -O2 -g -debug")
if(OSH_USE_OPENMP)
set(FLAGS "${FLAGS} -qopenmp")
endif()
else()
message(FATAL_ERROR "Unexpected compiler type ${CMAKE_CXX_COMPILER_ID}")
endif()
set(CMAKE_CXX_FLAGS "${FLAGS}")

configure_file(
  "${PROJECT_SOURCE_DIR}/omega_h_config_in.h"
  "${PROJECT_BINARY_DIR}/omega_h_config.h"
  )
include_directories(${PROJECT_BINARY_DIR})

add_library(omega_h omega_h.cpp)
if(OSH_USE_KOKKOS)
target_link_libraries(omega_h ${KOKKOS_LIB} dl)
endif()
if(OSH_USE_ZLIB)
target_link_libraries(omega_h ${ZLIB_LIBRARIES})
endif()

add_executable(unit_tests unit_tests.cpp)
target_link_libraries(unit_tests omega_h)
add_executable(perf_tests perf_tests.cpp)
target_link_libraries(perf_tests omega_h)
add_executable(mpi_tests mpi_tests.cpp)
target_link_libraries(mpi_tests omega_h)
add_executable(msh2osh msh2osh.cpp)
target_link_libraries(msh2osh omega_h)
add_executable(osh2vtk osh2vtk.cpp)
target_link_libraries(osh2vtk omega_h)
add_executable(oshdiff oshdiff.cpp)
target_link_libraries(oshdiff omega_h)
add_executable(ex1 ex1.cpp)
target_link_libraries(ex1 omega_h)
add_executable(ex2 ex2.cpp)
target_link_libraries(ex2 omega_h)
add_executable(ex3 ex3.cpp)
target_link_libraries(ex3 omega_h)
add_executable(flounder flounder.cpp)
target_link_libraries(flounder omega_h)
add_executable(corner_test corner_test.cpp)
target_link_libraries(corner_test omega_h)
add_executable(ring_test ring_test.cpp)
target_link_libraries(ring_test omega_h)
add_executable(coarsen_test coarsen_test.cpp)
target_link_libraries(coarsen_test omega_h)
add_executable(warp_test warp_test.cpp)
target_link_libraries(warp_test omega_h)

install(FILES
  "${PROJECT_BINARY_DIR}/omega_h_config.h"
  "${PROJECT_SOURCE_DIR}/omega_h.h"
  "${PROJECT_SOURCE_DIR}/omega_h.hpp"
  DESTINATION include)
install(TARGETS omega_h EXPORT omega_h-targets DESTINATION lib)
install(EXPORT omega_h-targets DESTINATION lib)
configure_file(
  "${PROJECT_SOURCE_DIR}/omega_h-config-in.cmake"
  "${PROJECT_BINARY_DIR}/omega_h-config.cmake"
  @ONLY
  )
install(FILES
  "${PROJECT_BINARY_DIR}/omega_h-config.cmake"
  DESTINATION lib/omega_h)
install(TARGETS
  msh2osh
  osh2vtk
  oshdiff
  DESTINATION bin)
