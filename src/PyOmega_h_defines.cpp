#include <Omega_h_defines.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_defines(py::module& module) {
  py::enum_<Omega_h_EntDim>(module, "EntDim", py::arithmetic())
      .value("VERT", OMEGA_H_VERT)
      .value("EDGE", OMEGA_H_EDGE)
      .value("FACE", OMEGA_H_FACE)
      .value("REGION", OMEGA_H_REGION)
      .export_values();
  py::enum_<Omega_h_Parting>(module, "Parting")
      .value("ELEM_BASED", OMEGA_H_ELEM_BASED)
      .value("GHOSTED", OMEGA_H_GHOSTED)
      .value("VERT_BASED", OMEGA_H_VERT_BASED)
      .export_values();
  py::enum_<Omega_h_Source>(
      module, "Source", "The type of source of a metric field")
      .value("CONSTANT", OMEGA_H_CONSTANT)
      .value("VARIATION", OMEGA_H_VARIATION)
      .value("DERIVATIVE", OMEGA_H_DERIVATIVE)
      .value("GIVEN", OMEGA_H_GIVEN)
      .value("IMPLIED", OMEGA_H_IMPLIED)
      //.value("PROXIMITY", OMEGA_H_PROXIMITY)
      .value("CURVATURE", OMEGA_H_CURVATURE)
      .export_values();
  py::enum_<Omega_h_Isotropy>(module, "Isotropy",
      "Whether and how to convert an anisotropic metric into an isotropic one")
      .value("ANISOTROPIC", OMEGA_H_ANISOTROPIC)
      .value("ISO_LENGTH", OMEGA_H_ISO_LENGTH)
      .value("ISO_SIZE", OMEGA_H_ISO_SIZE)
      .export_values();
  py::enum_<Omega_h_Scales>(module, "Scales",
      "Whether a metric source scales to satisfy element counts")
      .value("ABSOLUTE", OMEGA_H_ABSOLUTE)
      .value("SCALES", OMEGA_H_SCALES)
      .export_values();
  py::enum_<Omega_h_Family>(
      module, "Family", "Whether elements are simplices or hypercubes")
      .value("SIMPLEX", OMEGA_H_SIMPLEX)
      .value("HYPERCUBE", OMEGA_H_HYPERCUBE)
      .export_values();
}

}  // namespace Omega_h
