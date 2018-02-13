#include <Omega_h_array.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_array(py::module& module) {
  py::class_<Read<Real>>(module, "Read_real");
  py::class_<Reals, Read<Real>>(module, "Reals");
}

}  // namespace Omega_h
