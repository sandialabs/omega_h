#include <Omega_h_pybind11.hpp>

PYBIND11_MODULE(omega_h_pybind11, module) {
  Omega_h::pybind11_comm(module);
  Omega_h::pybind11_library(module);
}
