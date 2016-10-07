set(Omega_h_SOURCES
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
  derive.cpp
  project.cpp
  buffered_graphs.cpp
  overshoot.cpp
  collapse_vertex.cpp
  )

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

set(Gmodel_REQUIRED_VERSION 1.2.0)
set(Omega_h_USE_Gmodel_DEFAULT OFF)
bob_private_dep(Gmodel)

if(BUILD_TESTING)
  if(Omega_h_USE_MPI)
    get_filename_component(COMPILER_DIR ${CMAKE_CXX_COMPILER} PATH)
    find_program(MPIRUN mpirun PATHS ${COMPILER_DIR})
  endif()
  find_program(GMSH gmsh PATHS ${Gmsh_PREFIX})

  function(test_func TEST_NAME NUM_PROCS)
    if(MPIRUN)
      set(TEST_STR ${MPIRUN} -np ${NUM_PROCS} ${ARGN})
    else()
      if(NOT (${NUM_PROCS} EQUAL "1"))
        message(FATAL_ERROR "test_func with ${NUM_PROCS} procs but MPI disabled")
      endif()
      set(TEST_STR ${ARGN})
    endif()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_STR})
  endfunction(test_func)

  if(Omega_h_USE_Gmodel)
    function(gen_func MODEL_NAME MODEL_DIM)
      add_executable(gen_${MODEL_NAME} gen/${MODEL_NAME}.cpp)
      target_link_libraries(gen_${MODEL_NAME} PRIVATE gmodel)
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
    gen_func(twin_rotor 3)
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
  osh_add_exe(conserve_test)
  set(TEST_EXES ${TEST_EXES} conserve_test)
  test_func(run_conserve_test 1 ./conserve_test)
  osh_add_exe(perf_tests)
  osh_add_exe(ring_test)
  if(Omega_h_USE_Gmodel AND GMSH)
    set(TEST_EXES ${TEST_EXES} ring_test)
    test_func(run_ring_test 1 ./ring_test)
    set(TEST_MSHS ${TEST_MSHS} ring.msh)
  endif()
  osh_add_exe(aniso_test)
  test_func(run_aniso_test 1 ./aniso_test)
  osh_add_exe(imr_cases)
  osh_add_exe(test_normals)
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
      DEPENDS ${TEST_EXES} ${TEST_MSHS})
  osh_add_exe(test_degree)
endif()

set(SCOREC_REQUIRED_VERSION 2.0.0)
set(Omega_h_USE_SCOREC_DEFAULT OFF)
bob_private_dep(SCOREC)
if(Omega_h_USE_SCOREC)
  if(NOT Omega_h_USE_MPI)
    message(FATAL_ERROR "Building smb2osh requires enabling MPI")
  endif()
  osh_add_util(smb2osh)
  target_link_libraries(smb2osh PRIVATE SCOREC::core)
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
  Omega_h_mpi.h
  Omega_h_c.h
  Omega_h_kokkos.hpp
  Omega_h.hpp
  Omega_h_math.hpp
  Omega_h_functors.hpp
  Omega_h_r3d.hpp
  )

install(FILES ${Omega_h_HEADERS} DESTINATION include)

set(Meshb_REQUIRED_VERSION 7.12)
set(Omega_h_USE_Meshb_DEFAULT OFF)
bob_private_dep(Meshb)
if(Omega_h_USE_Meshb)
  osh_add_util(meshb2osh)
  target_link_libraries(meshb2osh PRIVATE Meshb::meshb)
endif()

bob_end_subdir()
