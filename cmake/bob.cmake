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
  message(STATUS "${PROJECT_NAME}_SEMVER = ${${PROJECT_NAME}_SEMVER}")
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

macro(bob_private_dep pkg_name)
  option(${PROJECT_NAME}_USE_${pkg_name} "Whether to use ${pkg_name}"
         ${${PROJECT_NAME}_USE_${pkg_name}_DEFAULT})
  message(STATUS "${PROJECT_NAME}_USE_${pkg_name}: ${${PROJECT_NAME}_USE_${pkg_name}}")
  if(${PROJECT_NAME}_USE_${pkg_name})
    set(${pkg_name}_PREFIX "${${pkg_name}_PREFIX_DEFAULT}"
        CACHE PATH "${pkg_name} install directory")
    if (${pkg_name}_PREFIX)
      message(STATUS "${pkg_name}_PREFIX ${${pkg_name}_PREFIX}")
      #if ${pkg_name}_PREFIX is set, don't find it anywhere else:
      find_package(${pkg_name} ${${pkg_name}_REQUIRED_VERSION}
                   REQUIRED PATHS ${${pkg_name}_PREFIX} NO_DEFAULT_PATH)
    else()
      #allow CMake to search other prefixes if ${pkg_name}_PREFIX is not set
      find_package(${pkg_name} ${${pkg_name}_REQUIRED_VERSION} REQUIRED)
    endif()
    if(${pkg_name}_CONFIG)
      message(STATUS "${pkg_name}_CONFIG: ${${pkg_name}_CONFIG}")
    endif()
    if(${pkg_name}_VERSION)
      message(STATUS "${pkg_name}_VERSION: ${${pkg_name}_VERSION}")
    endif()
  endif()
endmacro(bob_private_dep)

macro(bob_public_dep pkg_name)
  bob_private_dep(${pkg_name} "${version}")
  if(${PROJECT_NAME}_USE_${pkg_name})
    if (${pkg_name}_PREFIX)
      set(${PROJECT_NAME}_DEP_PREFIXES ${${PROJECT_NAME}_DEP_PREFIXES}
          ${${pkg_name}_PREFIX})
    endif()
    set(${PROJECT_NAME}_DEPS ${${PROJECT_NAME}_DEPS} ${pkg_name})
  endif()
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
  set(CONFIG_CONTENT
"set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION})
include(CMakeFindDependencyMacro)
# we will use find_dependency, but we don't want to force
# our users to have to specify where all of our dependencies
# were installed; that defeats the whole point of automatically
# importing dependencies.
# since the documentation for find_dependency() doesn't mention
# a PATHS argument, we'll temporarily add the prefixes to
# CMAKE_PREFIX_PATH.
set(${PROJECT_NAME}_DEPS \"${${PROJECT_NAME}_DEPS}\")
set(${PROJECT_NAME}_DEP_PREFIXES \"${${PROJECT_NAME}_DEP_PREFIXES}\")
set(${PROJECT_NAME}_BACKUP_PREFIX_PATH \"\${CMAKE_PREFIX_PATH}\")
set(CMAKE_PREFIX_PATH \"\${${PROJECT_NAME}_DEP_PREFIXES};\${CMAKE_PREFIX_PATH}\")
foreach(dep IN LISTS ${PROJECT_NAME}_DEPS)
  find_dependency(\${dep})
endforeach()
set(CMAKE_PREFIX_PATH \"\${${PROJECT_NAME}_BACKUP_PREFIX_PATH}\")
set(${PROJECT_NAME}_EXPORTED_TARGETS \"${${PROJECT_NAME}_EXPORTED_TARGETS}\")
foreach(tgt IN LISTS ${PROJECT_NAME}_EXPORTED_TARGETS)
  include(\${CMAKE_CURRENT_LIST_DIR}/\${tgt}-target.cmake)
endforeach()")
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
