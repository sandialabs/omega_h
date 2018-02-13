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
  py::class_<MetricSource>(module, "MetricSource", "Describes a single source metric field")
    .def(py::init<Omega_h_Source, Real, std::string const&,
                  Omega_h_Isotropy, Omega_h_Scales>(),
         "Build a MetricSource",
         py::arg("type"),
         py::arg("knob") = 1.0,
         py::arg("tag_name") = "",
         py::arg("isotropy") = OMEGA_H_ANISOTROPIC,
         py::arg("scales") = OMEGA_H_SCALES);
}

}  // namespace Omega_h


