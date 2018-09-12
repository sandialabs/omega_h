#ifndef OMEGA_H_PY_HPP
#define OMEGA_H_PY_HPP

#include <Omega_h_config.h>

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
void pybind11_defines(py::module& module);
void pybind11_array(py::module& module);
void pybind11_comm(py::module& module);
void pybind11_library(py::module& module);
void pybind11_mesh(py::module& module);
void pybind11_build(py::module& module);
void pybind11_adapt(py::module& module);
void pybind11_file(py::module& module);
void pybind11_class(py::module& module);
#ifdef OMEGA_H_USE_DOLFIN
void pybind11_dolfin(py::module& module);
#endif
}  // namespace Omega_h

#endif
