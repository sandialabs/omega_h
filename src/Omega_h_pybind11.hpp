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
void pybind11_comm(py::module& module);
void pybind11_library(py::module& module);
}  // namespace Omega_h

#endif
