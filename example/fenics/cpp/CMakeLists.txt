cmake_minimum_required(VERSION 3.0.1)

project(Omega_h-DOLFIN LANGUAGES C CXX)

find_package(DOLFIN REQUIRED)
# DOLFIN does conveniently provide us with its
# FindPETSc and FindSLEPC modules alongside
# its CMake configuration file
get_filename_component(DOLFIN_MODULES_DIR "${DOLFIN_CONFIG}" DIRECTORY)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${DOLFIN_MODULES_DIR}")

set(OMEGA_H_PREFIX "" CACHE PATH "Path to installed Omega_h")
find_package(Omega_h REQUIRED 9.9.0 PATHS ${OMEGA_H_PREFIX} NO_DEFAULT_PATH)

add_executable(omega_h-dolfin main.cpp)
target_link_libraries(omega_h-dolfin Omega_h::omega_h)
