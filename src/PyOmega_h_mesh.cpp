#include <Omega_h_mesh.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

#define OMEGA_H_DECL_TYPE(T, name)                                             \
  void (Mesh::*add_tag_##name)(Int, std::string const&, Int, Read<T>, bool) =  \
      &Mesh::add_tag<T>;
#define OMEGA_H_DEF_TYPE(T, name)                                              \
  .def("get_array_" #name, &Mesh::get_array<T>)                                \
      .def("add_tag_" #name, add_tag_##name,                                   \
          "Add " #name " tag array to the mesh",                               \
          py::arg("ent_dim") = OMEGA_H_VERT, py::arg("name"),                  \
          py::arg("ncomps") = 1, py::arg("array"),                             \
          py::arg("internal_do_not_use_ever") = false)

void pybind11_mesh(py::module& module) {
  OMEGA_H_DECL_TYPE(I8, int8)
  OMEGA_H_DECL_TYPE(I32, int32)
  OMEGA_H_DECL_TYPE(I64, int64)
  OMEGA_H_DECL_TYPE(Real, float64)
  void (Mesh::*set_parting)(Omega_h_Parting, Int, bool) = &Mesh::set_parting;
  void (Mesh::*balance)(bool) = &Mesh::balance;
  py::class_<Omega_h::Mesh>(module, "Mesh")
      .def("dim", &Omega_h::Mesh::dim)
      .def("nents", &Omega_h::Mesh::nents)
      .def("nglobal_ents", &Omega_h::Mesh::nglobal_ents)
      .def("set_parting", set_parting, py::arg("parting"), py::arg("nlayers"),
          py::arg("verbose") = false) OMEGA_H_DEF_TYPE(I8, int8)
          OMEGA_H_DEF_TYPE(I32, int32) OMEGA_H_DEF_TYPE(I64, int64)
              OMEGA_H_DEF_TYPE(Real, float64)
      .def("min_quality", &Omega_h::Mesh::min_quality)
      .def("max_length", &Omega_h::Mesh::max_length)
      .def("balance", balance, py::arg("predictive") = false);
  module.def(
      "new_empty_mesh", []() { return Mesh(pybind11_global_library.get()); });
}

#undef OMEGA_H_DECL_TYPE
#undef OMEGA_H_DEF_TYPE

}  // namespace Omega_h
