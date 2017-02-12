set(Omega_h_SOURCES
  Omega_h_c.cpp
  algebra.cpp
  Omega_h_control.cpp
  timer.cpp
  Omega_h_array.cpp
  int128.cpp
  repro.cpp
  sort.cpp
  scan.cpp
  indset.cpp
  map.cpp
  simplices.cpp
  form_uses.cpp
  find_unique.cpp
  invert_adj.cpp
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
  transfer_conserve.cpp
  Omega_h_compare.cpp
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
  derive.cpp
  project.cpp
  buffered_graphs.cpp
  Omega_h_coarsen_overshoot.cpp
  collapse_rail.cpp
  histogram.cpp
  eigen.cpp
  Omega_h_cmdline.cpp
  Omega_h_confined.cpp
  Omega_h_proximity.cpp
  Omega_h_motion_choices.cpp
  Omega_h_motion_transfer.cpp
  Omega_h_motion.cpp
  )

set(libMeshb_REQUIRED_VERSION 7.24)
set(Omega_h_USE_libMeshb_DEFAULT OFF)
bob_public_dep(libMeshb)
if(Omega_h_USE_libMeshb)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} meshb.cpp)
endif()

if(Omega_h_USE_EGADS)
  message(STATUS "Omega_h_USE_EGADS: ${Omega_h_USE_EGADS}")
  if(NOT EGADS_PREFIX)
    message(FATAL_ERROR "Omega_h_USE_EGADS=${Omega_h_USE_EGADS} but EGADS_PREFIX not set!")
  endif()
  find_path(EGADS_INCLUDE_DIR NAMES egads.h
      PATHS
      "${EGADS_PREFIX}/EngSketchPad/include"
      "${EGADS_PREFIX}/include"
      NO_DEFAULT_PATH)
  if(NOT EGADS_INCLUDE_DIR)
    message(FATAL_ERROR "could not find path to \"egads.h\"")
  endif()
  find_library(EGADS_LIBRARY NAMES egads
      PATHS
      "${EGADS_PREFIX}/EngSketchPad/lib"
      NO_DEFAULT_PATH)
  if(NOT EGADS_LIBRARY)
    message(FATAL_ERROR "could not find EGADS_LIBRARY")
  endif()
  message(STATUS "EGADS_INCLUDE_DIR: ${EGADS_INCLUDE_DIR}")
  message(STATUS "EGADS_LIBRARY: ${EGADS_LIBRARY}")
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_egads.cpp)
endif()

set(Omega_h_USE_SEACASExodus_DEFAULT OFF)
bob_public_dep(SEACASExodus)
if(Omega_h_USE_SEACASEXodus)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_exodus.cpp)
endif()

if(Omega_h_ONE_FILE)
  foreach(source IN LISTS Omega_h_SOURCES)
    set(ONE_FILE_CONTENTS "${ONE_FILE_CONTENTS}#include <${source}>\n")
  endforeach()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/omega_h.cpp "${ONE_FILE_CONTENTS}")
  add_library(omega_h omega_h.cpp)
else()
  add_library(omega_h ${Omega_h_SOURCES})
endif()

target_include_directories(omega_h
    PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
#for users to find all our headers after we're installed:
target_include_directories(omega_h INTERFACE $<INSTALL_INTERFACE:include>)
#find our generated omega_h_config.h, during build only
target_include_directories(omega_h
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

if(Omega_h_USE_libMeshb)
  target_link_libraries(omega_h PUBLIC libMeshb::Meshb.7)
endif()

if(Omega_h_USE_EGADS)
  target_include_directories(omega_h PUBLIC "${EGADS_INCLUDE_DIR}")
  target_link_libraries(omega_h PUBLIC "${EGADS_LIBRARY}")
endif()

if(Omega_h_USE_SEACASExodus)
  target_link_libraries(omega_h PUBLIC exodus)
endif()

bob_export_target(omega_h)

set(Omega_h_USE_Kokkos_DEFAULT OFF)
bob_public_dep(Kokkos)
if(Omega_h_USE_Kokkos)
  get_target_property(kokkoscore_LOCATION kokkoscore LOCATION)
  message(STATUS "kokkoscore LOCATION: ${kokkoscore_LOCATION}")
  get_target_property(kokkoscore_LINK kokkoscore INTERFACE_LINK_LIBRARIES)
  if(kokkoscore_LINK)
    message(STATUS "kokkoscore INTERFACE_LINK_LIBRARIES: ${kokkoscore_LINK}")
  endif()
  target_link_libraries(omega_h PUBLIC kokkoscore)
  target_include_directories(omega_h PUBLIC ${KokkosCore_INCLUDE_DIRS})
endif()

set(Omega_h_USE_ZLIB_DEFAULT ON)
bob_public_dep(ZLIB)
if(Omega_h_USE_ZLIB)
  message(STATUS "ZLIB_INCLUDE_DIRS: ${ZLIB_INCLUDE_DIRS}")
  # The ZLIB::ZLIB imported target is not present in CMake 3.0.0
  if (CMAKE_VERSION VERSION_LESS "3.1.0")
    target_link_libraries(omega_h PUBLIC ${ZLIB_LIBRARIES})
    target_include_directories(omega_h PUBLIC ${ZLIB_INCLUDE_DIRS})
  else()
    target_link_libraries(omega_h PUBLIC ZLIB::ZLIB)
  endif()
endif()

function(osh_add_exe EXE_NAME)
  add_executable(${EXE_NAME} ${EXE_NAME}.cpp)
  target_link_libraries(${EXE_NAME} PRIVATE omega_h)
endfunction(osh_add_exe)

macro(osh_add_util EXE_NAME)
  osh_add_exe(${EXE_NAME})
  bob_export_target(${EXE_NAME})
endmacro(osh_add_util)

osh_add_util(msh2osh)
osh_add_util(osh2vtk)
osh_add_util(oshdiff)
osh_add_util(vtkdiff)
osh_add_util(osh_box)
osh_add_util(osh_part)
osh_add_util(osh_scale)
if(Omega_h_USE_libMeshb)
  osh_add_util(meshb2osh)
  osh_add_util(osh2meshb)
  osh_add_util(ugawg_stats)
endif()

set(Gmodel_REQUIRED_VERSION 2.1.0)
set(Omega_h_USE_Gmodel_DEFAULT OFF)
bob_private_dep(Gmodel)

if(BUILD_TESTING)
  if(Omega_h_USE_MPI)
    get_filename_component(COMPILER_DIR ${CMAKE_CXX_COMPILER} PATH)
    find_program(MPIRUN mpirun PATHS ${COMPILER_DIR}
        DOC "Location of the mpirun program")
    if(NOT MPIRUN)
      message(STATUS "USE_MPI and BUILD_TESTING but MPIRUN not found!")
    endif()
  endif()

  function(test_func TEST_NAME NUM_PROCS)
    if(MPIRUN)
      set(TEST_STR ${MPIRUN} -np ${NUM_PROCS} ${ARGN})
    else()
      if(NOT (${NUM_PROCS} EQUAL "1"))
        message(STATUS "test ${TEST_NAME} ignored because MPIRUN not found!")
        return()
      endif()
      set(TEST_STR ${ARGN})
    endif()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_STR})
  endfunction(test_func)

  if(Omega_h_USE_Gmodel)
    function(gen_func MODEL_NAME)
      add_executable(gen_${MODEL_NAME} gen/${MODEL_NAME}.cpp)
      target_link_libraries(gen_${MODEL_NAME} PRIVATE gmodel)
      if(Omega_h_DATA)
        set(TEST_EXES ${TEST_EXES} gen_${MODEL_NAME})
        add_test(NAME run_gen_${MODEL_NAME} COMMAND ./gen_${MODEL_NAME})
        add_test(NAME gold_gen_${MODEL_NAME}
                 COMMAND diff
                 ./${MODEL_NAME}.geo
                 ${Omega_h_DATA}/${MODEL_NAME}.geo)
      endif()
    endfunction(gen_func)
    gen_func(ring)
    gen_func(ball_in_cube)
    gen_func(balls_in_box)
    gen_func(cylinder_thru_tube)
    gen_func(ball)
    gen_func(twin_rotor)
    gen_func(solder_balls)
  endif()

  osh_add_exe(unit_tests)
  set(TEST_EXES ${TEST_EXES} unit_tests)
  test_func(run_unit_tests 1 ./unit_tests)
  osh_add_exe(r3d_test)
  test_func(run_test_r3d 1 ./r3d_test)
  osh_add_exe(mpi_tests)
  set(TEST_EXES ${TEST_EXES} mpi_tests)
  if(Omega_h_USE_MPI)
    test_func(run_mpi_tests 4 ./mpi_tests)
  else()
    test_func(run_mpi_tests 1 ./mpi_tests)
  endif()
  osh_add_exe(corner_test)
  test_func(run_corner_test 1 ./corner_test)
  osh_add_exe(coarsen_test)
  test_func(run_coarsen_test 1 ./coarsen_test)
  osh_add_exe(warp_test)
  set(TEST_EXES ${TEST_EXES} warp_test)
  test_func(warp_test_serial 1 ./warp_test)
  if(Omega_h_USE_MPI)
    test_func(warp_test_parallel 2 ./warp_test)
  endif()
  osh_add_exe(2d_conserve_test)
  test_func(serial_2d_conserve 1 ./2d_conserve_test)
  if(Omega_h_USE_MPI)
    test_func(parallel_2d_conserve 2 ./2d_conserve_test)
  endif()
  osh_add_exe(3d_conserve_test)
  if(Omega_h_DATA)
    set(TEST_EXES ${TEST_EXES} 3d_conserve_test)
    test_func(serial_3d_conserve 1 ./3d_conserve_test
              ${Omega_h_DATA}/ball_in_cube.msh)
    if(Omega_h_USE_MPI)
      test_func(parallel_3d_conserve 2 ./3d_conserve_test
                ${Omega_h_DATA}/ball_in_cube.msh)
    endif()
  endif()
  osh_add_exe(perf_tests)
  osh_add_exe(ring_test)
  if(Omega_h_DATA)
    set(TEST_EXES ${TEST_EXES} ring_test)
    test_func(run_ring_test 1 ./ring_test ${Omega_h_DATA}/ring.msh)
  endif()
  osh_add_exe(aniso_test)
  test_func(run_aniso_test 1 ./aniso_test)
  osh_add_exe(imr_cases)
  if(Omega_h_DATA)
    osh_add_exe(curvature_test)
    set(TEST_EXES ${TEST_EXES} curvature_test)
    test_func(ring_curvature 1 ./curvature_test
              ${Omega_h_DATA} ring)
    test_func(ball_curvature 1 ./curvature_test
              ${Omega_h_DATA} ball)
    if(Omega_h_USE_MPI)
      test_func(parallel_ring_curvature 2 ./curvature_test
                ${Omega_h_DATA} ring)
      test_func(parallel_ball_curvature 2 ./curvature_test
                ${Omega_h_DATA} ball)
    endif()
  endif()
  osh_add_exe(bend_test)
  osh_add_exe(test_degree)
  if(Omega_h_USE_libMeshb)
    osh_add_exe(ugawg_linear)
    if(Omega_h_USE_EGADS)
      osh_add_exe(ugawg_cylinder)
    endif()
  endif()
  osh_add_exe(slope_test)
  if(Omega_h_USE_EGADS)
    osh_add_exe(egads_test)
  endif()
  osh_add_exe(advect2d_test)
  if(Omega_h_DATA)
    set(TEST_EXES ${TEST_EXES} advect2d_test)
    test_func(advect2d_serial 1 ./advect2d_test ${Omega_h_DATA}/advect2d.osh)
    test_func(advect2d_parallel 2 ./advect2d_test ${Omega_h_DATA}/advect2d.osh)
  endif()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
      DEPENDS ${TEST_EXES})
endif()

set(DEF_VARS
    VERSION_MAJOR
    VERSION_MINOR
    VERSION_PATCH
    USE_MPI
    USE_Kokkos
    USE_OpenMP
    USE_CUDA
    USE_ZLIB
    USE_libMeshb
    USE_EGADS
    USE_SEACASExodus
    CHECK_BOUNDS
    PROTECT
    )
foreach(def_var IN LISTS DEF_VARS)
  string(TOUPPER "Omega_h_${def_var}" uppercase_var)
  set(${uppercase_var} ${Omega_h_${def_var}})
endforeach()
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Omega_h_config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h")

set(Omega_h_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h"
  Omega_h_macros.h
  Omega_h_mpi.h
  Omega_h_c.h
  Omega_h_kokkos.hpp
  Omega_h_defines.hpp
  Omega_h_array.hpp
  Omega_h.hpp
  Omega_h_math.hpp
  Omega_h_functors.hpp
  )

install(FILES ${Omega_h_HEADERS} DESTINATION include)

bob_end_subdir()
