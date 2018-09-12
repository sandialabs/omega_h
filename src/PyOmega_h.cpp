#include <PyOmega_h.hpp>

PYBIND11_MODULE(PyOmega_h, module) {
  module.doc() = "Omega_h: simplex mesh adaptation";
  Omega_h::pybind11_defines(module);
  Omega_h::pybind11_array(module);
  Omega_h::pybind11_comm(module);
  Omega_h::pybind11_library(module);
  Omega_h::pybind11_mesh(module);
  Omega_h::pybind11_build(module);
  Omega_h::pybind11_adapt(module);
  Omega_h::pybind11_file(module);
  Omega_h::pybind11_class(module);
#ifdef OMEGA_H_USE_DOLFIN
  Omega_h::pybind11_dolfin(module);
#endif
}
