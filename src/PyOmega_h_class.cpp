#include <Omega_h_class.hpp>
#include <Omega_h_mesh.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_class(py::module& module) {
  module.def("classify_by_angles", &Omega_h::classify_by_angles,
      "Classify a mesh, sharp-angle creases become corners/curves",
      py::arg("mesh"), py::arg("sharp_angle"));
}

}  // namespace Omega_h
