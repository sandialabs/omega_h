#include <Omega_h_build.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_build(py::module& module) {
  module.def("build_box", &Omega_h::build_box, "Build a rectangular mesh",
      // if we give this a default argument, its lifetime is that
      // of a global variable, which exceeds the Library lifetime!
      py::arg("comm") /*= pybind11_global_library->world()*/,
      py::arg("family") = OMEGA_H_SIMPLEX, py::arg("x") = 1.0,
      py::arg("y") = 1.0, py::arg("z") = 1.0, py::arg("nx") = 0,
      py::arg("ny") = 0, py::arg("nz") = 0, py::arg("symmetric") = false);
}

}  // namespace Omega_h
