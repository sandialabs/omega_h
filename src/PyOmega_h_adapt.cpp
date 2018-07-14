#include <Omega_h_adapt.hpp>
#include <Omega_h_mesh.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_adapt(py::module& module) {
  py::enum_<Verbosity>(
      module, "Verbosity", "Level of output verbosity in adaptation")
      .value("SILENT", SILENT)
      .value("EACH_ADAPT", EACH_ADAPT)
      .value("EACH_REBUILD", EACH_REBUILD)
      .value("EXTRA_STATS", EXTRA_STATS)
      .export_values();
  py::class_<AdaptOpts>(
      module, "AdaptOpts", "Options controlling adaptation behavior")
      .def(py::init<Mesh*>())
      .def_readwrite("verbosity", &AdaptOpts::verbosity)
      .def_readwrite("min_quality_allowed", &AdaptOpts::min_quality_allowed);
  py::class_<MetricSource>(
      module, "MetricSource", "Describes a single source metric field")
      .def(py::init<Omega_h_Source, Real, std::string const&, Omega_h_Isotropy,
               Omega_h_Scales>(),
          "Build a MetricSource", py::arg("type"), py::arg("knob") = 1.0,
          py::arg("tag_name") = "", py::arg("isotropy") = OMEGA_H_ANISOTROPIC,
          py::arg("scales") = OMEGA_H_SCALES);
  py::class_<MetricInput>(
      module, "MetricInput", "Describes all inputs that form a metric field")
      .def(py::init<>())
      .def_readwrite("verbose", &MetricInput::verbose)
      .def_readwrite("should_limit_lengths", &MetricInput::should_limit_lengths)
      .def_readwrite("max_length", &MetricInput::max_length)
      .def_readwrite("min_length", &MetricInput::min_length)
      .def_readwrite(
          "should_limit_gradation", &MetricInput::should_limit_gradation)
      .def_readwrite("max_gradation_rate", &MetricInput::max_gradation_rate)
      .def_readwrite("gradation_convergence_tolerance",
          &MetricInput::gradation_convergence_tolerance)
      .def_readwrite("should_limit_element_count",
          &MetricInput::should_limit_element_count)
      .def_readwrite("max_element_count", &MetricInput::max_element_count)
      .def_readwrite("min_element_count", &MetricInput::min_element_count)
      .def_readwrite("element_count_over_relaxation",
          &MetricInput::element_count_over_relaxation)
      .def_readwrite("nsmoothing_steps", &MetricInput::nsmoothing_steps)
      .def("add_source", &MetricInput::add_source);
  module.def("grade_fix_adapt", &grade_fix_adapt,
      "Apply gradation control, possibly fix quality, then adapt",
      py::arg("mesh"), py::arg("opts"), py::arg("target_metric"),
      py::arg("verbose") = true);
  module.def("add_implied_metric_tag", &add_implied_metric_tag);
  module.def("generate_target_metric_tag", &generate_target_metric_tag);
  module.def("approach_metric", &approach_metric, py::arg("mesh"),
      py::arg("opts"), py::arg("min_step") = 1e-4);
  module.def("adapt", &adapt);
}

}  // namespace Omega_h
