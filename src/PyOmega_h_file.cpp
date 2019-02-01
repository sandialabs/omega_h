#include <Omega_h_file.hpp>
#include <Omega_h_filesystem.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_file(py::module& module) {
  py::class_<Omega_h::filesystem::path>(module, "path")
      .def(py::init<char const*>());
  Mesh (*gmsh_read_file)(filesystem::path const&, CommPtr) = &gmsh::read;
  void (*gmsh_write_file)(filesystem::path const&, Mesh*) = &gmsh::write;
  void (*vtk_write_vtu_dim)(std::string const&, Mesh*, Int, bool) =
      &vtk::write_vtu;
  void (*vtk_write_vtu)(std::string const&, Mesh*, bool) = &vtk::write_vtu;
  void (*vtk_write_parallel_dim)(std::string const&, Mesh*, Int, bool) =
      &vtk::write_parallel;
  void (*vtk_write_parallel)(std::string const&, Mesh*, bool) =
      &vtk::write_parallel;
  module.def("gmsh_read_file", gmsh_read_file, "Read a Gmsh file");
  module.def("gmsh_write_file", gmsh_write_file, "Write a Gmsh file");
  module.def("vtk_write_vtu", vtk_write_vtu, "Write a mesh as a .vtu file",
      py::arg("path"), py::arg("mesh"), py::arg("compress") = true);
  module.def("vtk_write_vtu_dim", vtk_write_vtu_dim,
      "Write entities of one dimension as a .vtu file", py::arg("path"),
      py::arg("mesh"), py::arg("cell_dim"), py::arg("compress") = true);
  module.def("vtk_write_parallel", vtk_write_parallel,
      "Write a mesh as a directory of parallel VTK files", py::arg("path"),
      py::arg("mesh"), py::arg("compress") = true);
  module.def("vtk_write_parallel_dim", vtk_write_parallel_dim,
      "Write entities of one dimension as a directory of parallel VTK files",
      py::arg("path"), py::arg("mesh"), py::arg("cell_dim"),
      py::arg("compress") = true);
}

}  // namespace Omega_h
