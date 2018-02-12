#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <pybind11/pybind11.h>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#include <Omega_h_library.hpp>

namespace py = pybind11;

PYBIND11_MODULE(omega_h_pybind11, module) {
  py::class_<Omega_h::Library>(module, "Library")
    .def(py::init<>())
    .def("world", &Omega_h::Library::world);
}
