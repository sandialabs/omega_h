#include <Omega_h_array.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

void pybind11_array(py::module& module) {
  py::class_<Write<Real>>(module, "Write_float64")
    .def("size", &Write<Real>::size)
    ;
  py::class_<Read<Real>>(module, "Read_float64")
    .def(py::init<Write<Real>>())
    .def("size", &Read<Real>::size)
    ;
  py::class_<Reals, Read<Real>>(module, "Reals")
    .def(py::init<LO, Real, std::string const&>(),
        py::arg("size"),
        py::arg("value"),
        py::arg("name") = "")
    ;
  py::class_<HostRead<Real>>(module, "HostRead_float64", py::buffer_protocol())
    .def(py::init<Read<Real>>())
    .def_buffer([](HostRead<Real>& a) -> py::buffer_info {
        return py::buffer_info(
  /* note: although this breaks const-correctness, I think
     the benefits may outweigh the drawbacks */
            const_cast<Real*>(a.data()),
            sizeof(Real),
            py::format_descriptor<Real>::format(),
            1,
            { a.size() },
            { sizeof(Real) }
            );
        })
     .def("get", &HostRead<Real>::get)
     ;
  py::class_<HostWrite<Real>>(module, "HostWrite_float64", py::buffer_protocol())
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
     .def("write", &HostWrite<Real>::write)
     ;
  Write<Real> (*deep_copy_float64)(Read<Real> a)
    = &deep_copy;
  module.def("deep_copy_float64", deep_copy_float64);
}

}  // namespace Omega_h
