#include <Omega_h_mesh.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_mesh(py::module& module) {
  void (Mesh::*add_tag_real)(Int, std::string const&, Int, Read<Real>, bool)
    = &Mesh::add_tag<Real>;
  py::class_<Omega_h::Mesh>(module, "Mesh")
      .def("dim", &Omega_h::Mesh::dim)
      .def("nents", &Omega_h::Mesh::nents)
      .def("nglobal_ents", &Omega_h::Mesh::nglobal_ents)
      .def("get_array_real", &Omega_h::Mesh::get_array<Omega_h::Real>)
      .def("add_tag_real", add_tag_real,
          "Add a tag array to the mesh",
          py::arg("ent_dim") = OMEGA_H_VERT,
          py::arg("name"),
          py::arg("ncomps") = 1,
          py::arg("array"),
          py::arg("internal_do_not_use_ever") = false)
      ;
}

}  // namespace Omega_h
