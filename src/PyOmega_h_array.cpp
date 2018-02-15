#include <Omega_h_array.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_array(py::module& module) {
  py::class_<Write<Real>>(module, "Write_real")
    .def("size", &Write<Real>::size)
    ;
  py::class_<Read<Real>>(module, "Read_real")
    .def(py::init<Write<Real>>())
    .def("size", &Read<Real>::size)
    ;
  py::class_<Reals, Read<Real>>(module, "Reals")
    .def(py::init<LO, Real, std::string const&>(),
        py::arg("size"),
        py::arg("value"),
        py::arg("name") = "")
    ;
  py::class_<HostWrite<Real>>(module, "HostWrite_real", py::buffer_protocol())
    .def(py::init<LO, std::string const&>(),
        py::arg("size"),
        py::arg("name") = "")
    .def_buffer([](HostWrite<Real>& a) -> py::buffer_info {
        return py::buffer_info(
            a.data(),
            sizeof(Real),
            py::format_descriptor<Real>::format(),
            1,
            { a.size() },
            { sizeof(Real) }
            );
        })
     .def("set", &HostWrite<Real>::set)
     .def("get", &HostWrite<Real>::get)
     ;
  Write<Real> (*deep_copy_real)(Read<Real> a)
    = &deep_copy;
  module.def("deep_copy_real", deep_copy_real);
}

}  // namespace Omega_h
