function(bob_always_full_rpath)
  # CMake RPATH "always full" configuration, see:
  # https://cmake.org/Wiki/CMake_RPATH_handling#Always_full_RPATH
  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH False PARENT_SCOPE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH False PARENT_SCOPE)
  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
       "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" PARENT_SCOPE)
  endif()
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True PARENT_SCOPE)
endfunction(bob_always_full_rpath)

macro(bob_begin_package)
  message(STATUS "CMAKE_VERSION: ${CMAKE_VERSION}")
  if (${PROJECT_NAME}_VERSION)
    message(STATUS "${PROJECT_NAME}_VERSION: ${${PROJECT_NAME}_VERSION}")
  endif()
  #try to force BUILD_TESTING to be OFF by default
  set(BUILD_TESTING OFF CACHE BOOL "Build and run tests")
  include(CTest)
  enable_testing()
  option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
  #If not building shared libs, then prefer static
  #dependency libs
  if(NOT BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".dylib")
  endif()
  bob_always_full_rpath()
  message(STATUS "BUILD_TESTING: ${BUILD_TESTING}")
  message(STATUS "BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")
  message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
endmacro(bob_begin_package)

function(bob_form_semver)
  execute_process(COMMAND git describe --exact-match HEAD
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE NOT_TAG
      OUTPUT_VARIABLE TAG_NAME
      ERROR_VARIABLE TAG_ERROR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
  if(NOT_TAG)
    if(${PROJECT_NAME}_VERSION)
      set(SEMVER ${${PROJECT_NAME}_VERSION})
      execute_process(COMMAND git log -1 --format=%h
          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
          RESULT_VARIABLE NO_SHA1
          OUTPUT_VARIABLE SHORT_SHA1
          ERROR_VARIABLE SHA1_ERROR
          OUTPUT_STRIP_TRAILING_WHITESPACE
          )
      if(NO_SHA1)
        message(WARNING "bob_form_semver no Git hash !\n" ${SHA1_ERROR})
      else()
        set(SEMVER "${SEMVER}-sha.${SHORT_SHA1}")
      endif()
    else()
      message(FATAL_ERROR "bob_form_semver needs either ${PROJECT_NAME}_VERSION or a Git tag\n" ${TAG_ERROR})
    endif()
  else()
    if(TAG_NAME MATCHES "^v([0-9]+[.])?([0-9]+[.])?([0-9]+)$")
      string(SUBSTRING "${TAG_NAME}" 1 -1 SEMVER)
      if(${PROJECT_NAME}_VERSION AND (NOT (SEMVER VERSION_EQUAL ${PROJECT_NAME}_VERSION)))
        message(FATAL_ERROR "bob_form_semver: tag is ${TAG_NAME} but ${PROJECT_NAME}_VERSION=${${PROJECT_NAME}_VERSION} !")
      endif()
    else()
      if(${PROJECT_NAME}_VERSION)
        set(SEMVER "${${PROJECT_NAME}_VERSION}-tag.${TAG_NAME}")
      else()
        message(FATAL_ERROR "bob_form_semver needs either ${PROJECT_NAME}_VERSION or a Git tag of the form v1.2.3")
      endif()
    endif()
  endif()
  if(${PROJECT_NAME}_KEY_BOOLS)
    set(SEMVER "${SEMVER}+")
    foreach(KEY_BOOL IN LISTS ${PROJECT_NAME}_KEY_BOOLS)
      if(${KEY_BOOL})
        set(SEMVER "${SEMVER}1")
      else()
        set(SEMVER "${SEMVER}0")
      endif()
    endforeach()
  endif()
  set(${PROJECT_NAME}_SEMVER "${SEMVER}" PARENT_SCOPE)
  message(STATUS "${PROJECT_NAME}_SEMVER = ${SEMVER}")
endfunction(bob_form_semver)

function(bob_begin_cxx_flags)
  if(CMAKE_BUILD_TYPE)
    message(FATAL_ERROR "can't set CMAKE_BUILD_TYPE and use bob_*_cxx_flags")
  endif()
  option(${PROJECT_NAME}_CXX_OPTIMIZE "Compile C++ with optimization" ON)
  option(${PROJECT_NAME}_CXX_SYMBOLS "Compile C++ with debug symbols" ON)
  option(${PROJECT_NAME}_CXX_WARNINGS "Compile C++ with warnings" ON)
  set(FLAGS "")
  if(${PROJECT_NAME}_CXX_OPTIMIZE)
    set(FLAGS "${FLAGS} -O3")
  else()
    set(FLAGS "${FLAGS} -O0")
  endif()
  if(${PROJECT_NAME}_CXX_SYMBOLS)
    set(FLAGS "${FLAGS} -g")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (${PROJECT_NAME}_CXX_WARNINGS)
      set(FLAGS "${FLAGS} -Werror -Weverything")
      set(FLAGS "${FLAGS} -Wno-padded")
      set(FLAGS "${FLAGS} -Wno-float-equal")
      set(FLAGS "${FLAGS} -Wno-weak-template-vtables")
    endif()
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    if (${PROJECT_NAME}_CXX_WARNINGS)
      set(FLAGS "${FLAGS} -Werror -Wall -Wextra")
      set(FLAGS "${FLAGS} -Wlogical-op -Wold-style-cast")
      set(FLAGS "${FLAGS} -Wdouble-promotion -Wshadow -Wformat=2")
      if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "6.0")
        set(FLAGS "${FLAGS} -Wduplicated-cond -Wnull-dereference")
      endif()
      if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.0")
        set(FLAGS "${FLAGS} -Wduplicated-branches -Wrestrict")
      endif()
    endif()
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  else()
    message(WARNING "Unexpected compiler type ${CMAKE_CXX_COMPILER_ID}")
  endif()
  set(CMAKE_CXX_FLAGS "${FLAGS}" PARENT_SCOPE)
endfunction(bob_begin_cxx_flags)

function(bob_cxx11_flags)
  set(FLAGS "${CMAKE_CXX_FLAGS}")
  set(FLAGS "${FLAGS} --std=c++11")
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if (${PROJECT_NAME}_CXX_WARNINGS)
      set(FLAGS "${FLAGS} -Wno-c++98-compat-pedantic -Wno-c++98-compat")
    endif()
  endif()
  set(CMAKE_CXX_FLAGS "${FLAGS}" PARENT_SCOPE)
endfunction(bob_cxx11_flags)

function(bob_end_cxx_flags)
  set(${PROJECT_NAME}_CXX_FLAGS "" CACHE STRING "Override all C++ compiler flags")
  set(${PROJECT_NAME}_EXTRA_CXX_FLAGS "" CACHE STRING "Extra C++ compiler flags")
  set(FLAGS "${CMAKE_CXX_FLAGS}")
  if(${PROJECT_NAME}_CXX_FLAGS)
    set(FLAGS "${${PROJECT_NAME}_CXX_FLAGS}")
  else()
    set(FLAGS "${FLAGS} ${${PROJECT_NAME}_EXTRA_CXX_FLAGS}")
  endif()
  message(STATUS "CMAKE_CXX_FLAGS: ${FLAGS}")
  set(CMAKE_CXX_FLAGS "${FLAGS}" PARENT_SCOPE)
endfunction(bob_end_cxx_flags)

macro(bob_add_dependency)
  set(options PUBLIC PRIVATE)
  set(oneValueArgs NAME)
  set(multiValueArgs COMPONENTS)
  cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  if (NOT ARG_NAME)
    message(FATAL_ERROR "bob_add_dependency: no NAME argument given")
  endif()
  if (ARG_PUBLIC AND ARG_PRIVATE)
    message(FATAL_ERROR "bob_add_dependency: can't specify both PUBLIC and PRIVATE")
  endif()
  if (ARG_COMPONENTS)
    set(ARG_COMPONENTS COMPONENTS ${ARG_COMPONENTS})
  endif()
  option(${PROJECT_NAME}_USE_${ARG_NAME} "Whether to use ${ARG_NAME}"
         ${${PROJECT_NAME}_USE_${ARG_NAME}_DEFAULT})
  message(STATUS "${PROJECT_NAME}_USE_${ARG_NAME}: ${${PROJECT_NAME}_USE_${ARG_NAME}}")
  if(${PROJECT_NAME}_USE_${ARG_NAME})
    set(${ARG_NAME}_PREFIX "${${ARG_NAME}_PREFIX_DEFAULT}"
        CACHE PATH "${ARG_NAME} install directory")
    if (${ARG_NAME}_PREFIX)
      message(STATUS "${ARG_NAME}_PREFIX ${${ARG_NAME}_PREFIX}")
      #if ${ARG_NAME}_PREFIX is set, don't find it anywhere else:
      set(ARG_PREFIX PATHS "${${ARG_NAME}_PREFIX}" NO_DEFAULT_PATH)
    else()
      #allow CMake to search other prefixes if ${ARG_NAME}_PREFIX is not set
      set(ARG_PREFIX)
    endif()
    set(${ARG_NAME}_find_package_args
        "${${ARG_NAME}_REQUIRED_VERSION}"
        ${ARG_COMPONENTS}
        ${ARG_PREFIX})
    find_package(${ARG_NAME} ${${ARG_NAME}_find_package_args} REQUIRED)
    if(${ARG_NAME}_CONFIG)
      message(STATUS "${ARG_NAME}_CONFIG: ${${ARG_NAME}_CONFIG}")
    endif()
    if(${ARG_NAME}_VERSION)
      message(STATUS "${ARG_NAME}_VERSION: ${${ARG_NAME}_VERSION}")
    endif()
  endif()
  if (ARG_PUBLIC AND ${PROJECT_NAME}_USE_${ARG_NAME})
    set(${PROJECT_NAME}_DEPS ${${PROJECT_NAME}_DEPS} ${ARG_NAME})
  endif()
endmacro(bob_add_dependency)

macro(bob_private_dep pkg_name)
  bob_add_dependency(PRIVATE NAME "${pkg_name}")
endmacro(bob_private_dep)

macro(bob_public_dep pkg_name)
  bob_add_dependency(PUBLIC NAME "${pkg_name}")
endmacro(bob_public_dep)

function(bob_target_includes lib_name)
  #find local headers even with #include <>
  target_include_directories(${lib_name}
      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
  #find generated configuration headers
  target_include_directories(${lib_name}
      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
endfunction(bob_target_includes)

function(bob_library_includes lib_name)
  bob_target_includes("${lib_name}")
  #ensure downstream users include installed headers
  target_include_directories(${lib_name} INTERFACE $<INSTALL_INTERFACE:include>)
endfunction(bob_library_includes)

function(bob_export_target tgt_name)
  install(TARGETS ${tgt_name} EXPORT ${tgt_name}-target
      RUNTIME DESTINATION bin
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib)
  install(EXPORT ${tgt_name}-target NAMESPACE ${PROJECT_NAME}::
          DESTINATION lib/cmake/${PROJECT_NAME})
  set(${PROJECT_NAME}_EXPORTED_TARGETS
      ${${PROJECT_NAME}_EXPORTED_TARGETS} ${tgt_name} PARENT_SCOPE)
endfunction(bob_export_target)

macro(bob_end_subdir)
  set(${PROJECT_NAME}_EXPORTED_TARGETS
      ${${PROJECT_NAME}_EXPORTED_TARGETS} PARENT_SCOPE)
  set(${PROJECT_NAME}_DEPS ${${PROJECT_NAME}_DEPS} PARENT_SCOPE)
  set(${PROJECT_NAME}_DEP_PREFIXES ${${PROJECT_NAME}_DEP_PREFIXES} PARENT_SCOPE)
endmacro(bob_end_subdir)

function(bob_config_header HEADER_PATH)
  get_filename_component(HEADER_NAME "${HEADER_PATH}" NAME)
  string(REPLACE "." "_" INCLUDE_GUARD "${HEADER_NAME}")
  string(TOUPPER "${INCLUDE_GUARD}" INCLUDE_GUARD)
  set(HEADER_CONTENT
"#ifndef ${INCLUDE_GUARD}
#define ${INCLUDE_GUARD}
")
  if (${PROJECT_NAME}_KEY_BOOLS)
    foreach(KEY_BOOL IN LISTS ${PROJECT_NAME}_KEY_BOOLS)
      if (${KEY_BOOL})
        string(TOUPPER "${KEY_BOOL}" MACRO_NAME)
        set(HEADER_CONTENT
"${HEADER_CONTENT}
#define ${MACRO_NAME}")
      endif()
    endforeach()
  endif()
  if (${PROJECT_NAME}_KEY_INTS)
    foreach(KEY_INT IN LISTS ${PROJECT_NAME}_KEY_INTS)
      string(TOUPPER "${KEY_INT}" MACRO_NAME)
      set(HEADER_CONTENT
"${HEADER_CONTENT}
#define ${MACRO_NAME} ${${KEY_INT}}")
    endforeach()
  endif()
  if (${PROJECT_NAME}_KEY_STRINGS)
    foreach(KEY_STRING IN LISTS ${PROJECT_NAME}_KEY_STRINGS)
      string(TOUPPER "${KEY_STRING}" MACRO_NAME)
      set(HEADER_CONTENT
"${HEADER_CONTENT}
#define ${MACRO_NAME} \"${${KEY_STRING}}\"")
    endforeach()
  endif()
  set(HEADER_CONTENT
"${HEADER_CONTENT}

#endif
")
  file(WRITE "${HEADER_PATH}" "${HEADER_CONTENT}")
endfunction()

function(bob_end_package)
  include(CMakePackageConfigHelpers)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR lib)
  set(LATEST_FIND_DEPENDENCY
"#The definition of this macro is really inconvenient prior to CMake
#commit ab358d6a859d8b7e257ed1e06ca000e097a32ef6
#we'll just copy the latest code into our Config.cmake file
macro(latest_find_dependency dep)
  if (NOT \${dep}_FOUND)
    set(cmake_fd_quiet_arg)
    if(\${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY)
      set(cmake_fd_quiet_arg QUIET)
    endif()
    set(cmake_fd_required_arg)
    if(\${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
      set(cmake_fd_required_arg REQUIRED)
    endif()

    get_property(cmake_fd_alreadyTransitive GLOBAL PROPERTY
      _CMAKE_\${dep}_TRANSITIVE_DEPENDENCY
    )

    find_package(\${dep} \${ARGN}
      \${cmake_fd_quiet_arg}
      \${cmake_fd_required_arg}
    )

    if(NOT DEFINED cmake_fd_alreadyTransitive OR cmake_fd_alreadyTransitive)
      set_property(GLOBAL PROPERTY _CMAKE_\${dep}_TRANSITIVE_DEPENDENCY TRUE)
    endif()

    if (NOT \${dep}_FOUND)
      set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE \"\${CMAKE_FIND_PACKAGE_NAME} could not be found because dependency \${dep} could not be found.\")
      set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND False)
      return()
    endif()
    set(cmake_fd_required_arg)
    set(cmake_fd_quiet_arg)
    set(cmake_fd_exact_arg)
  endif()
endmacro(latest_find_dependency)"
       )
  set(FIND_DEPS_CONTENT)
  foreach(dep IN LISTS ${PROJECT_NAME}_DEPS)
    string(REPLACE ";" " " FIND_DEP_ARGS "${${dep}_find_package_args}")
    set(FIND_DEPS_CONTENT
"${FIND_DEPS_CONTENT}
latest_find_dependency(${dep} ${FIND_DEP_ARGS})"
       )
  endforeach()
  set(CONFIG_CONTENT
"set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION})
${LATEST_FIND_DEPENDENCY}
${FIND_DEPS_CONTENT}
set(${PROJECT_NAME}_EXPORTED_TARGETS \"${${PROJECT_NAME}_EXPORTED_TARGETS}\")
foreach(tgt IN LISTS ${PROJECT_NAME}_EXPORTED_TARGETS)
  include(\${CMAKE_CURRENT_LIST_DIR}/\${tgt}-target.cmake)
endforeach()"
  )
  foreach(TYPE IN ITEMS "BOOL" "INT" "STRING")
    if (${PROJECT_NAME}_KEY_${TYPE}S)
      foreach(KEY_${TYPE} IN LISTS ${PROJECT_NAME}_KEY_${TYPE}S)
        set(CONFIG_CONTENT
"${CONFIG_CONTENT}
set(${KEY_${TYPE}} \"${${KEY_${TYPE}}}\")")
      endforeach()
    endif()
  endforeach()
  set(CONFIG_CONTENT
"${CONFIG_CONTENT}
")
  install(FILES
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    DESTINATION lib/cmake/${PROJECT_NAME})
  if(PROJECT_VERSION)
    file(WRITE
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        "${CONFIG_CONTENT}")
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)
    install(FILES
      "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
      DESTINATION lib/cmake/${PROJECT_NAME})
  endif()
endfunction(bob_end_package)
