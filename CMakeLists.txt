cmake_minimum_required(VERSION 2.6)
project(omega_h CXX)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
set(CMAKE_CXX_FLAGS "\
--std=c++11 \
-O2 \
-g \
-fno-omit-frame-pointer \
-ferror-limit=1 \
-Werror \
-Weverything \
-Wno-c++98-compat \
-Wno-c++98-compat-pedantic \
-Wno-padded \
-Wno-float-equal \
")
else ()
  message(FATAL_ERROR "Unexpected compiler type ${CMAKE_CXX_COMPILER_ID}")
endif()

#configure_file(to_config.hpp config.hpp)
#include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library(omega_h omega_h.cpp)
add_executable(unit_tests unit_tests.cpp)
target_link_libraries(unit_tests omega_h)
# add_executable(perf_tests perf_tests.cpp)
# target_link_libraries(perf_tests omega_h)
