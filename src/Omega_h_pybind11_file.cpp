#include <Omega_h_file.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_file(py::module& module) {
  Mesh (*gmsh_read_file)(std::string const&, CommPtr) = &gmsh::read;
  void (*gmsh_write_file)(std::string const&, Mesh*) = &gmsh::write;
  module.def("gmsh_read_file", gmsh_read_file, "Read a Gmsh file");
  module.def("gmsh_write_file", gmsh_write_file, "Read a Gmsh file");
}

}
