#include <Omega_h_dolfin.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_dolfin(py::module& module) {
  void (*mesh_from_dolfin)(Mesh*, dolfin::Mesh const&)
    = &from_dolfin;
  void (*function_from_dolfin)(Mesh*, dolfin::Function const&, std::string const&)
    = &from_dolfin;
  module.def("mesh_to_dolfin", &to_dolfin,
      "Convert an Omega_h mesh to a DOLFIN mesh");
  module.def("mesh_from_dolfin", &mesh_from_dolfin,
      "Convert a DOLFIN mesh to an Omega_h mesh");
  module.def("function_from_dolfin", &function_from_dolfin,
      "Convert a DOLFIN Function to Omega_h tag(s)");
}

}
