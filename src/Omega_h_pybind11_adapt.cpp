#include <Omega_h_adapt.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_pybind11.hpp>

namespace Omega_h {

void pybind11_adapt(py::module& module) {
  py::enum_<Verbosity>(module, "Verbosity", "Level of output verbosity in adaptation")
    .value("SILENT", SILENT)
    .value("EACH_ADAPT", EACH_ADAPT)
    .value("EACH_REBUILD", EACH_REBUILD)
    .value("EXTRA_STATS", EXTRA_STATS)
    .export_values();
  py::class_<AdaptOpts>(module, "AdaptOpts", "Options controlling adaptation behavior")
    .def(py::init<Mesh*>())
    .def_readwrite("verbosity", &AdaptOpts::verbosity);
}

}  // namespace Omega_h


