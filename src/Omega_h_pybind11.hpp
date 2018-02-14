#ifndef OMEGA_H_PYBIND11_HPP
#define OMEGA_H_PYBIND11_HPP

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <pybind11/pybind11.h>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

namespace py = pybind11;

namespace Omega_h {
class Library;
extern std::unique_ptr<Library> pybind11_global_library;
void pybind11_c(py::module& module);
void pybind11_array(py::module& module);
void pybind11_comm(py::module& module);
void pybind11_library(py::module& module);
void pybind11_mesh(py::module& module);
void pybind11_build(py::module& module);
void pybind11_adapt(py::module& module);
void pybind11_file(py::module& module);
}  // namespace Omega_h

#endif
