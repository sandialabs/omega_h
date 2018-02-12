#include <Omega_h_mesh.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_mesh(py::module& module) {
  py::class_<Omega_h::Mesh>(module, "Mesh")
      .def("dim", &Omega_h::Mesh::dim)
      .def("nents", &Omega_h::Mesh::nents)
      .def("nglobal_ents", &Omega_h::Mesh::nglobal_ents)
      ;
}

}  // namespace Omega_h
