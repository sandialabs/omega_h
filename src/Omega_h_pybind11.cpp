#include <Omega_h_pybind11.hpp>

PYBIND11_MODULE(omega_h_pybind11, module) {
  module.doc() = "Omega_h: simplex mesh adaptation";
  Omega_h::pybind11_c(module);
  Omega_h::pybind11_array(module);
  Omega_h::pybind11_comm(module);
  Omega_h::pybind11_library(module);
  Omega_h::pybind11_mesh(module);
  Omega_h::pybind11_build(module);
  Omega_h::pybind11_adapt(module);
  Omega_h::pybind11_file(module);
#ifdef OMEGA_H_USE_DOLFIN
  Omega_h::pybind11_dolfin(module);
#endif
}
