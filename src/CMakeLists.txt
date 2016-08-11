configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/omega_h_config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/omega_h_config.h")

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
  )
endif()

#for users to find all our headers after we're installed:
target_include_directories(omega_h INTERFACE $<INSTALL_INTERFACE:include>)

#find our generated omega_h_config.h, during build only
target_include_directories(omega_h
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

macro(osh_find_package PKG_NAME VERSION)
  set(${PKG_NAME}_PREFIX "" CACHE PATH "${PKG_NAME} install directory")
  if ("${${PKG_NAME}_PREFIX}" STREQUAL "")
    #allow CMake to search other prefixes if Kokkos_PREFIX is not set
    find_package(${PKG_NAME} ${VERSION} REQUIRED)
  else()
    #but if Kokkos_PREFIX is set, don't find it anywhere else:
    find_package(${PKG_NAME} ${VERSION} REQUIRED PATHS ${${PKG_NAME}_PREFIX} NO_DEFAULT_PATH)
  endif()
  message(STATUS "${PKG_NAME}_CONFIG: ${${PKG_NAME}_CONFIG}")
endmacro(osh_find_package)

if(OSH_USE_KOKKOS)
  osh_find_package(Kokkos "")
  get_target_property(kokkoscore_LOCATION kokkoscore LOCATION)
  message(STATUS "kokkoscore LOCATION: ${kokkoscore_LOCATION}")
  get_target_property(kokkoscore_LINK kokkoscore INTERFACE_LINK_LIBRARIES)
  message(STATUS "kokkoscore INTERFACE_LINK_LIBRARIES: ${kokkoscore_LINK}")
  target_link_libraries(omega_h PUBLIC kokkoscore)
  target_include_directories(omega_h PUBLIC ${KokkosCore_INCLUDE_DIRS})
  set(DEP_LIBS ${DEP_LIBS} ${KOKKOS_LIB})
  set(DEP_INCLUDE_DIRS ${DEP_INCLUDE_DIRS} ${KOKKOS_INCLUDE})
endif()

if(OSH_USE_ZLIB)
  if(NOT "${ZLIB_PREFIX}" STREQUAL "")
    set(ZLIB_ROOT "${ZLIB_PREFIX}")
  endif()
  find_package(ZLIB REQUIRED)
  message(STATUS "ZLIB_INCLUDE_DIRS: ${ZLIB_INCLUDE_DIRS}")
  target_link_libraries(omega_h PRIVATE ${ZLIB_LIBRARIES})
  target_include_directories(omega_h PRIVATE ${ZLIB_INCLUDE_DIRS})
endif()

function(exe_func EXE_NAME)
  add_executable(${EXE_NAME} ${EXE_NAME}.cpp)
  target_link_libraries(${EXE_NAME} PRIVATE omega_h)
endfunction(exe_func)

function(util_func EXE_NAME)
  exe_func(${EXE_NAME})
  install(TARGETS ${EXE_NAME} EXPORT ${EXE_NAME}-targets DESTINATION bin)
  install(EXPORT ${EXE_NAME}-targets NAMESPACE Omega_h:: DESTINATION lib/omega_h)
  set(OSH_TARGETS ${OSH_TARGETS} ${EXE_NAME} PARENT_SCOPE)
endfunction(util_func)

util_func(msh2osh)
util_func(osh2vtk)
util_func(oshdiff)
util_func(vtkdiff)

if(BUILD_TESTING)
  if(OSH_USE_MPI)
    get_filename_component(COMPILER_DIR ${CMAKE_CXX_COMPILER} PATH)
    find_program(MPIRUN mpirun PATHS ${COMPILER_DIR})
  endif()
  find_program(GMSH gmsh PATHS ${Gmsh_PREFIX})
  find_package(Gmodel 1.2.0 QUIET PATHS ${Gmodel_PREFIX})

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

  if(Gmodel_FOUND)
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

  exe_func(unit_tests)
  set(TEST_EXES ${TEST_EXES} unit_tests)
  test_func(run_unit_tests 1 ./unit_tests)
  exe_func(mpi_tests)
  set(TEST_EXES ${TEST_EXES} mpi_tests)
  if(OSH_USE_MPI)
    test_func(run_mpi_tests 4 ./mpi_tests)
  else()
    test_func(run_mpi_tests 1 ./mpi_tests)
  endif()
  exe_func(warp_test)
  set(TEST_EXES ${TEST_EXES} warp_test)
  test_func(warp_test_serial 1 ./warp_test)
  if(OSH_USE_MPI)
    test_func(warp_test_parallel 2 ./warp_test)
  endif()
  exe_func(perf_tests)
  exe_func(flounder)
  exe_func(corner_test)
  exe_func(coarsen_test)
  exe_func(aniso_test)
  exe_func(ring_test)
  exe_func(imr_cases)
  exe_func(test_normals)
  if(Gmodel_FOUND AND GMSH)
    set(TEST_EXES ${TEST_EXES} ring_test)
    test_func(run_ring_test 1 ./ring_test)
    set(TEST_MSHS ${TEST_MSHS} ring.msh)
  endif()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
      DEPENDS ${TEST_EXES} ${TEST_MSHS})
  if(OSH_USE_KOKKOS)
    exe_func(test_view)
  endif()
  exe_func(test_r3d)
  test_func(run_test_r3d 1 ./test_r3d)
  exe_func(test_hessian)
  exe_func(test_degree)
endif()

if(OSH_BUILD_SMB2OSH)
  if(NOT OSH_USE_MPI)
    message(FATAL_ERROR "Building smb2osh requires enabling MPI")
  endif()
  osh_find_package(SCOREC 2.0.0)
  util_func(smb2osh)
  target_link_libraries(smb2osh PRIVATE SCOREC::core)
  install(TARGETS smb2osh DESTINATION bin)
endif()

set(OSH_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/omega_h_config.h"
  omega_h_mpi.h
  omega_h_c.h
  omega_h_kokkos.hpp
  omega_h.hpp
  omega_h_math.hpp
  omega_h_functors.hpp
  omega_h_r3d.hpp
  )

install(FILES ${OSH_HEADERS} DESTINATION include)
install(TARGETS omega_h EXPORT omega_h-targets DESTINATION lib)
install(EXPORT omega_h-targets DESTINATION lib/omega_h)

set(OSH_TARGETS ${OSH_TARGETS} omega_h PARENT_SCOPE)
