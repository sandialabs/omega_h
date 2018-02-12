#include <Omega_h_build.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_build(py::module& module) {
  module.def("build_box", &Omega_h::build_box, "Build a rectangular mesh",
      py::arg("comm"),
      py::arg("family") = OMEGA_H_SIMPLEX,
      py::arg("x") = 1.0,
      py::arg("y") = 1.0,
      py::arg("z") = 1.0,
      py::arg("nx") = 0,
      py::arg("ny") = 0,
      py::arg("nz") = 0);
}

}  // namespace Omega_h

