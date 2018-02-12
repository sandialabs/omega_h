#include <Omega_h_c.h>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_c(py::module& module) {
  py::enum_<Omega_h_Family>(module, "Family", "Whether elements are simplices or hypercubes")
    .value("SIMPLEX", OMEGA_H_SIMPLEX)
    .value("HYPERCUBE", OMEGA_H_HYPERCUBE)
    .export_values();
}

}  // namespace Omega_h


