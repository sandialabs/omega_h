#include <Omega_h_pybind11.hpp>
#include <Omega_h_comm.hpp>

namespace Omega_h {

void pybind11_comm(py::module& module) {
  py::class_<Omega_h::Comm>(module, "Comm")
    .def("rank", &Omega_h::Comm::rank)
    .def("size", &Omega_h::Comm::size)
    ;
}

}
