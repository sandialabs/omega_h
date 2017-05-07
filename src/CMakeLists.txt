set(Omega_h_SOURCES
  Omega_h_c.cpp
  Omega_h_control.cpp
  Omega_h_timer.cpp
  Omega_h_array.cpp
  Omega_h_array_ops.cpp
  Omega_h_vector.cpp
  Omega_h_matrix.cpp
  Omega_h_map.cpp
  Omega_h_scan.cpp
  Omega_h_sort.cpp
  Omega_h_int128.cpp
  Omega_h_indset.cpp
  Omega_h_simplex.cpp
  Omega_h_adj.cpp
  Omega_h_tag.cpp
  Omega_h_mesh.cpp
  Omega_h_bbox.cpp
  Omega_h_hilbert.cpp
  Omega_h_build.cpp
  Omega_h_graph.cpp
  Omega_h_file.cpp
  Omega_h_base64.cpp
  Omega_h_vtk.cpp
  Omega_h_simplify.cpp
  Omega_h_box.cpp
  Omega_h_surface.cpp
  Omega_h_mark.cpp
  Omega_h_class.cpp
  Omega_h_reorder.cpp
  Omega_h_shape.cpp
  Omega_h_quality.cpp
  Omega_h_gmsh.cpp
  Omega_h_comm.cpp
  Omega_h_remotes.cpp
  Omega_h_dist.cpp
  Omega_h_linpart.cpp
  Omega_h_owners.cpp
  Omega_h_migrate.cpp
  Omega_h_bcast.cpp
  Omega_h_unmap_mesh.cpp
  Omega_h_ghost.cpp
  Omega_h_inertia.cpp
  Omega_h_bipart.cpp
  Omega_h_metric.cpp
  Omega_h_refine_qualities.cpp
  Omega_h_refine_topology.cpp
  Omega_h_modify.cpp
  Omega_h_refine.cpp
  Omega_h_transfer.cpp
  Omega_h_conserve.cpp
  Omega_h_compare.cpp
  Omega_h_check_collapse.cpp
  Omega_h_coarsen_qualities.cpp
  Omega_h_coarsen_topology.cpp
  Omega_h_coarsen.cpp
  Omega_h_approach.cpp
  Omega_h_laplace.cpp
  Omega_h_adapt.cpp
  Omega_h_swap.cpp
  Omega_h_swap2d_qualities.cpp
  Omega_h_swap2d_topology.cpp
  Omega_h_swap2d.cpp
  Omega_h_swap3d_qualities.cpp
  Omega_h_swap3d_topology.cpp
  Omega_h_swap3d.cpp
  Omega_h_regress.cpp
  Omega_h_consistent.cpp
  Omega_h_align.cpp
  Omega_h_xml.cpp
  Omega_h_recover.cpp
  Omega_h_coarsen_overshoot.cpp
  Omega_h_collapse_rail.cpp
  Omega_h_histogram.cpp
  Omega_h_eigen.cpp
  Omega_h_cmdline.cpp
  Omega_h_confined.cpp
  Omega_h_proximity.cpp
  Omega_h_motion_choices.cpp
  Omega_h_motion_transfer.cpp
  Omega_h_motion.cpp
  Omega_h_metric_input.cpp
  Omega_h_assoc.cpp
  )

if(Omega_h_USE_libMeshb)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_meshb.cpp)
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

if(Omega_h_USE_SEACASExodus)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_exodus.cpp)
endif()

if(USE_TEUCHOS)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_teuchos.cpp)
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

if(Omega_h_USE_Kokkos)
  target_link_libraries(omega_h PUBLIC kokkoscore)
  target_include_directories(omega_h PUBLIC ${KokkosCore_INCLUDE_DIRS})
endif()

if(Omega_h_USE_libMeshb)
  target_link_libraries(omega_h PUBLIC libMeshb::Meshb.7)
endif()

if(Omega_h_USE_EGADS)
  target_include_directories(omega_h PUBLIC "${EGADS_INCLUDE_DIR}")
  target_link_libraries(omega_h PUBLIC "${EGADS_LIBRARY}")
endif()

if(Omega_h_USE_SEACASExodus)
  target_link_libraries(omega_h PUBLIC exodus)
  target_include_directories(omega_h
      PRIVATE ${SEACASExodus_TPL_INCLUDE_DIRS})
  target_link_libraries(omega_h
      PUBLIC ${SEACASExodus_TPL_LIBRARIES})
endif()

if(Omega_h_USE_TeuchosParameterList)
  target_link_libraries(omega_h PUBLIC teuchosparameterlist)
  target_include_directories(omega_h
      PUBLIC ${TeuchosParameterList_TPL_INCLUDE_DIRS})
  target_link_libraries(omega_h
      PUBLIC ${TeuchosParameterList_TPL_LIBRARIES})
endif()

if(Omega_h_USE_TeuchosComm)
  target_link_libraries(omega_h PUBLIC teuchoscomm)
  target_include_directories(omega_h
      PUBLIC ${TeuchosComm_TPL_INCLUDE_DIRS})
  target_link_libraries(omega_h
      PUBLIC ${TeuchosComm_TPL_LIBRARIES})
endif()

bob_export_target(omega_h)

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
if(Omega_h_USE_SEACASExodus)
  osh_add_util(exo2osh)
  osh_add_util(osh2exo)
endif()
if(USE_TEUCHOS)
  osh_add_util(osh_adapt)
  osh_add_util(osh_copy_params)
endif()

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
  osh_add_exe(curvature_test)
  if(Omega_h_DATA)
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
    if(Omega_h_USE_MPI)
      test_func(advect2d_parallel 2 ./advect2d_test ${Omega_h_DATA}/advect2d.osh)
    endif()
  endif()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
      DEPENDS ${TEST_EXES})
  osh_add_exe(1d_test)
  test_func(serial_1d_test 1 ./1d_test)
  if(Omega_h_USE_MPI)
    test_func(parallel_1d_test 2 ./1d_test)
  endif()
endif()

bob_config_header("${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h")
set(Omega_h_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h"
  Omega_h_macros.h
  Omega_h_mpi.h
  Omega_h_c.h
  Omega_h_kokkos.hpp
  Omega_h_defines.hpp
  Omega_h_array.hpp
  Omega_h.hpp
  Omega_h_vector.hpp
  Omega_h_matrix.hpp
  Omega_h_functors.hpp
  Omega_h_array_ops.hpp
  Omega_h_graph.hpp
  Omega_h_map.hpp
  Omega_h_tag.hpp
  Omega_h_int128.hpp
  Omega_h_comm.hpp
  Omega_h_remotes.hpp
  Omega_h_dist.hpp
  Omega_h_adj.hpp
  Omega_h_library.hpp
  Omega_h_mesh.hpp
  Omega_h_file.hpp
  Omega_h_adapt.hpp
  Omega_h_build.hpp
  Omega_h_class.hpp
  Omega_h_scalar.hpp
  Omega_h_align.hpp
  Omega_h_few.hpp
  Omega_h_mark.hpp
  Omega_h_loop.hpp
  Omega_h_timer.hpp
  Omega_h_eigen.hpp
  Omega_h_lie.hpp
  Omega_h_recover.hpp
  Omega_h_metric.hpp
  Omega_h_assoc.hpp
  Omega_h_cmdline.hpp
  Omega_h_compare.hpp
  Omega_h_simplex.hpp
  Omega_h_qr.hpp
  Omega_h_shape.hpp
  Omega_h_quality.hpp
  )

if(USE_TEUCHOS)
  set(Omega_h_HEADERS ${Omega_h_HEADERS} Omega_h_teuchos.hpp)
endif()

install(FILES ${Omega_h_HEADERS} DESTINATION include)

bob_end_subdir()
