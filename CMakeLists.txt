cmake_minimum_required(VERSION 2.8.11)
project(Omega_h CXX)

#trying to force BUILD_TESTING to default to OFF
set(BUILD_TESTING OFF CACHE STRING "Build and run tests")

include(CTest)
enable_testing()

set(OSH_MAJOR "1")
set(OSH_MINOR "1")
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
set(SCOREC_PREFIX "" CACHE STRING "SCOREC/core install directory")
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(OSH_SANITIZE_ADDRESS "Use -fsanitize=address" OFF)
set(Kokkos_PREFIX "" CACHE STRING "Kokkos installation directory")
option(OSH_PROTECT "Catch OS signals and print stack" OFF)
set(Gmodel_PREFIX "" CACHE STRING "Gmodel installation directory")
set(Gmsh_PREFIX "" CACHE STRING "Gmsh installation directory")
option(OSH_ONE_FILE "Compile omega_h as a single file" ON)

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
    set(FLAGS "${FLAGS} --expt-extended-lambda")
    set(FLAGS "${FLAGS} -arch=sm_${OSH_CUDA_ARCH}")
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

if(OSH_ONE_FILE)
  add_library(omega_h omega_h.cpp)
else()
  add_library(omega_h
  algebra.cpp
  control.cpp
  protect.cpp
  timer.cpp
  array.cpp
  int128.cpp
  repro.cpp
  sort.cpp
  scan.cpp
  indset.cpp
  map.cpp
  simplices.cpp
  form_uses.cpp
  find_unique.cpp
  invert.cpp
  reflect_down.cpp
  transit.cpp
  tag.cpp
  mesh.cpp
  bbox.cpp
  hilbert.cpp
  construct.cpp
  graph.cpp
  star.cpp
  file.cpp
  base64.cpp
  vtk.cpp
  simplify.cpp
  box.cpp
  surface.cpp
  mark.cpp
  classify.cpp
  reorder.cpp
  size.cpp
  access.cpp
  quality.cpp
  gmsh.cpp
  comm.cpp
  remotes.cpp
  dist.cpp
  linpart.cpp
  owners.cpp
  migrate.cpp
  bcast.cpp
  unmap_mesh.cpp
  ghost.cpp
  inertia.cpp
  bipart.cpp
  metric.cpp
  refine_qualities.cpp
  refine_topology.cpp
  modify.cpp
  refine.cpp
  transfer.cpp
  compare.cpp
  check_collapse.cpp
  coarsen_qualities.cpp
  coarsen_topology.cpp
  coarsen.cpp
  approach.cpp
  laplace.cpp
  adapt.cpp
  swap.cpp
  swap2d_qualities.cpp
  swap2d_topology.cpp
  swap2d.cpp
  swap3d_qualities.cpp
  swap3d_topology.cpp
  swap3d.cpp
  regress.cpp
  consistent.cpp
  align.cpp
  xml.cpp
  )
endif()

target_link_libraries(omega_h ${DEP_LIBS})

add_executable(msh2osh msh2osh.cpp)
target_link_libraries(msh2osh omega_h)
add_executable(osh2vtk osh2vtk.cpp)
target_link_libraries(osh2vtk omega_h)
add_executable(oshdiff oshdiff.cpp)
target_link_libraries(oshdiff omega_h)
add_executable(vtkdiff vtkdiff.cpp)
target_link_libraries(vtkdiff omega_h)

if(BUILD_TESTING)
  if(OSH_USE_MPI)
    get_filename_component(COMPILER_DIR ${CMAKE_CXX_COMPILER} PATH)
    find_program(MPIRUN mpirun PATHS ${COMPILER_DIR})
  endif()
  find_program(GMSH gmsh PATHS ${Gmsh_PREFIX})
  find_package(Gmodel 1.0.1 QUIET PATHS ${Gmodel_PREFIX})

  function(test_func TEST_NAME NUM_PROCS)
    if(MPIRUN)
      set(TEST_STR ${MPIRUN} -np ${NUM_PROCS} ${ARGN})
    else()
      if(NOT (${NUM_PROCS} EQUAL "1"))
        message(FATAL_ERROR "test_func with ${NUM_PROCS} procs but MPI disabled")
      endif()
      set(TEST_STR ${MPIRUN} -np ${NUM_PROCS} ${ARGN})
    endif()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_STR})
  endfunction(test_func)

  if(Gmodel_FOUND)
    function(gen_func MODEL_NAME MODEL_DIM)
      add_executable(gen_${MODEL_NAME} gen/${MODEL_NAME}.cpp)
      target_include_directories(gen_${MODEL_NAME} PRIVATE
                                 ${Gmodel_INCLUDE_DIRS})
      target_link_libraries(gen_${MODEL_NAME} gmodel)
      add_custom_command(OUTPUT ${MODEL_NAME}.geo
                         COMMAND ./gen_${MODEL_NAME}
                         DEPENDS gen_${MODEL_NAME}
                         COMMENT "generating ${MODEL_NAME}.geo")
      if(GMSH)
        if("${MODEL_DIM}" EQUAL "2")
          set(GMSH_CMD ${GMSH} -2 -o ${MODEL_NAME}.msh ${MODEL_NAME}.geo)
        else()
          set(GMSH_CMD ${GMSH} -3 -optimize_netgen -o ${MODEL_NAME}.msh ${MODEL_NAME}.geo)
        endif()
        add_custom_command(OUTPUT ${MODEL_NAME}.msh
                           COMMAND ${GMSH_CMD}
                           DEPENDS ${MODEL_NAME}.geo
                           COMMENT "generating ${MODEL_NAME}.msh")
      endif()
    endfunction(gen_func)
    gen_func(ring 2)
    gen_func(ball_in_cube 3)
    gen_func(balls_in_box 3)
    gen_func(cylinder_thru_tube 3)
    gen_func(ball 3)
  endif()

  add_executable(unit_tests unit_tests.cpp)
  target_link_libraries(unit_tests omega_h)
  set(TEST_EXES ${TEST_EXES} unit_tests)
  test_func(run_unit_tests 1 ./unit_tests)
  add_executable(mpi_tests mpi_tests.cpp)
  target_link_libraries(mpi_tests omega_h)
  set(TEST_EXES ${TEST_EXES} mpi_tests)
  if(OSH_USE_MPI)
    test_func(run_mpi_tests 4 ./mpi_tests)
  else()
    test_func(run_mpi_tests 1 ./mpi_tests)
  endif()
  add_executable(warp_test warp_test.cpp)
  target_link_libraries(warp_test omega_h)
  set(TEST_EXES ${TEST_EXES} warp_test)
  test_func(warp_test_serial 1 ./warp_test)
  test_func(warp_test_parallel 2 ./warp_test)
  add_executable(perf_tests perf_tests.cpp)
  target_link_libraries(perf_tests omega_h)
  add_executable(flounder flounder.cpp)
  target_link_libraries(flounder omega_h)
  add_executable(corner_test corner_test.cpp)
  target_link_libraries(corner_test omega_h)
  add_executable(coarsen_test coarsen_test.cpp)
  target_link_libraries(coarsen_test omega_h)
  add_executable(aniso_test aniso_test.cpp)
  target_link_libraries(aniso_test omega_h)
  add_executable(ring_test ring_test.cpp)
  target_link_libraries(ring_test omega_h)
  add_executable(imr_cases imr_cases.cpp)
  target_link_libraries(imr_cases omega_h)
  add_executable(test_normals test_normals.cpp)
  target_link_libraries(test_normals omega_h)
  if(Gmodel_FOUND AND GMSH)
    set(TEST_EXES ${TEST_EXES} ring_test)
    test_func(run_ring_test 1 ./ring_test)
    set(TEST_MSHS ${TEST_MSHS} ring.msh)
  endif()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
      DEPENDS ${TEST_EXES} ${TEST_MSHS})
endif()

if(OSH_BUILD_SMB2OSH)
  find_path(SCOREC_INCLUDE apf.h PATHS ${SCOREC_PREFIX}/include)
  find_library(MDS_LIB NAMES mds PATHS ${SCOREC_PREFIX}/lib)
  find_library(APF_LIB NAMES apf PATHS ${SCOREC_PREFIX}/lib)
  find_library(MTH_LIB NAMES mth PATHS ${SCOREC_PREFIX}/lib)
  find_library(LION_LIB NAMES lion PATHS ${SCOREC_PREFIX}/lib)
  find_library(GMI_LIB NAMES gmi PATHS ${SCOREC_PREFIX}/lib)
  find_library(PCU_LIB NAMES pcu PATHS ${SCOREC_PREFIX}/lib)
  add_executable(smb2osh smb2osh.cpp)
  target_include_directories(smb2osh PRIVATE ${SCOREC_INCLUDE})
  target_link_libraries(smb2osh
      omega_h
      ${APF_LIB}
      ${MDS_LIB}
      ${GMI_LIB}
      ${PCU_LIB})
  if (MTH_LIB AND LION_LIB)
    target_link_libraries(smb2osh
        ${MTH_LIB}
        ${LION_LIB})
  endif()
  install(TARGETS smb2osh DESTINATION bin)
endif()

set(OSH_HEADERS
  "${PROJECT_BINARY_DIR}/omega_h_config.h"
  "${PROJECT_SOURCE_DIR}/omega_h_mpi.h"
  "${PROJECT_SOURCE_DIR}/omega_h_c.h"
  "${PROJECT_SOURCE_DIR}/omega_h_kokkos.hpp"
  "${PROJECT_SOURCE_DIR}/omega_h.hpp"
  )

install(FILES ${OSH_HEADERS} DESTINATION include)
install(TARGETS omega_h EXPORT omega_h-targets DESTINATION lib)
install(EXPORT omega_h-targets DESTINATION lib/omega_h)
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
install(TARGETS
  msh2osh
  osh2vtk
  oshdiff
  vtkdiff
  DESTINATION bin)
