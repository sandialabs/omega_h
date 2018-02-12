#include <Omega_h_pybind11.hpp>
#include <Omega_h_library.hpp>

PYBIND11_MODULE(omega_h_pybind11, module) {
  py::class_<Omega_h::Library>(module, "Library")
    .def(py::init<>())
    .def("world", &Omega_h::Library::world);
}
