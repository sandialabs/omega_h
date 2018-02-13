#include <Omega_h_array.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_array(py::module& module) {
  py::class_<Read<Real>>(module, "Read_Real");
//py::class_<Omega_h::Reals>(module, "Reals");
}

}  // namespace Omega_h

