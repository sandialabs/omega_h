cmake_minimum_required(VERSION 2.8.11)
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
option(OSH_BUILD_SMB2OSH "Enable building the smb2osh program" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(OSH_SANITIZE_ADDRESS "Use -fsanitize=address" OFF)
set(Kokkos_PREFIX "" CACHE STRING "Kokkos installation directory")
option(OSH_PROTECT "Catch OS signals and print stack" OFF)

if(OSH_USE_KOKKOS)
  find_path(KOKKOS_INCLUDE Kokkos_Core.hpp PATHS "${Kokkos_PREFIX}/include")
  find_library(KOKKOS_LIB NAMES kokkos kokkoscore PATHS "${Kokkos_PREFIX}/lib")
  set(DEP_LIBS ${DEP_LIBS} ${KOKKOS_LIB})
  set(DEP_INCLUDE_DIRS ${DEP_INCLUDE_DIRS} ${KOKKOS_INCLUDE})
endif()

if(OSH_USE_ZLIB)
  find_package(ZLIB REQUIRED)
  set(DEP_LIBS ${DEP_LIBS} ${ZLIB_LIBRARIES})
  set(DEP_INCLUDE_DIRS ${DEP_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
endif()

include_directories(${DEP_INCLUDE_DIRS})

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
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(FLAGS "${FLAGS} -fno-omit-frame-pointer -ferror-limit=1")
  set(FLAGS "${FLAGS} -Werror -Weverything -Wno-c++98-compat")
  set(FLAGS "${FLAGS} -Wno-c++98-compat-pedantic -Wno-padded")
  set(FLAGS "${FLAGS} -Wno-float-equal")
  set(FLAGS "${FLAGS} -Wno-weak-template-vtables")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  if(OSH_USE_CUDA)
    set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
    set(FLAGS "${FLAGS} -expt-extended-lambda -arch=sm_${OSH_CUDA_ARCH}")
  elseif(OSH_USE_OPENMP)
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer -fmax-errors=1")
    set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
    set(FLAGS "${FLAGS} -fopenmp")
  else()
    set(FLAGS "${FLAGS} -fno-omit-frame-pointer -fmax-errors=1")
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
set(CMAKE_CXX_FLAGS "${FLAGS}")

configure_file(
  "${PROJECT_SOURCE_DIR}/omega_h_config.h.in"
  "${PROJECT_BINARY_DIR}/omega_h_config.h")
include_directories(${PROJECT_BINARY_DIR})

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

add_library(omega_h omega_h.cpp)
target_link_libraries(omega_h ${DEP_LIBS})

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
add_executable(vtkdiff vtkdiff.cpp)
target_link_libraries(vtkdiff omega_h)
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
add_executable(imr_cases imr_cases.cpp)
target_link_libraries(imr_cases omega_h)
add_executable(aniso_test aniso_test.cpp)
target_link_libraries(aniso_test omega_h)

if(OSH_BUILD_SMB2OSH)
find_path(APF_INCLUDE apf.h)
find_library(MDS_LIB NAMES mds)
find_library(APF_LIB NAMES apf)
find_library(MTH_LIB NAMES mth)
find_library(LION_LIB NAMES lion)
find_library(GMI_LIB NAMES gmi)
find_library(PCU_LIB NAMES pcu)
add_executable(smb2osh smb2osh.cpp)
target_include_directories(smb2osh PRIVATE ${APF_INCLUDE})
target_link_libraries(smb2osh
    omega_h
    ${APF_LIB}
    ${MTH_LIB}
    ${LION_LIB}
    ${MDS_LIB}
    ${GMI_LIB}
    ${PCU_LIB})
install(TARGETS smb2osh DESTINATION bin)
endif()

install(FILES
  "${PROJECT_BINARY_DIR}/omega_h_config.h"
  "${PROJECT_SOURCE_DIR}/omega_h_mpi.h"
  "${PROJECT_SOURCE_DIR}/omega_h_c.h"
  "${PROJECT_SOURCE_DIR}/omega_h_kokkos.hpp"
  "${PROJECT_SOURCE_DIR}/omega_h.hpp"
  DESTINATION include)
install(TARGETS omega_h EXPORT omega_h-targets DESTINATION lib)
install(EXPORT omega_h-targets DESTINATION lib)
configure_file(
  "${PROJECT_SOURCE_DIR}/omega_h-config-in.cmake"
  "${PROJECT_BINARY_DIR}/omega_h-config.cmake"
  @ONLY)
configure_file(
  "${PROJECT_SOURCE_DIR}/omega_h-config-version-in.cmake"
  "${PROJECT_BINARY_DIR}/omega_h-config-version.cmake"
  @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}/omega_h-config.cmake"
  "${PROJECT_BINARY_DIR}/omega_h-config-version.cmake"
  DESTINATION lib/omega_h)
install(TARGETS
  msh2osh
  osh2vtk
  oshdiff
  vtkdiff
  DESTINATION bin)

configure_file(
  "${PROJECT_SOURCE_DIR}/run_tests.sh"
  "${PROJECT_BINARY_DIR}/run_tests.sh"
  COPYONLY)
