#include <Omega_h_pybind11.hpp>
#include <Omega_h_library.hpp>

namespace Omega_h {

void pybind11_library(py::module& module) {
  py::class_<Omega_h::Library>(module, "Library")
    .def(py::init<>())
    .def("world", &Omega_h::Library::world);
}

}
