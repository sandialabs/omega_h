#include <Omega_h_array.hpp>
#include <PyOmega_h.hpp>

namespace Omega_h {

template <class Scalar, class Wrapper>
static void pybind11_array_type(py::module& module,
    std::string const& py_scalar, std::string const& py_wrapper) {
  auto write_name = std::string("Write_") + py_scalar;
  auto read_name = std::string("Read_") + py_scalar;
  auto hostread_name = std::string("HostRead_") + py_scalar;
  auto hostwrite_name = std::string("HostWrite_") + py_scalar;
  auto deepcopy_name = std::string("deep_copy_") + py_scalar;
  py::class_<Write<Scalar>>(module, write_name.c_str())
      .def("size", &Write<Scalar>::size);
  py::class_<Read<Scalar>>(module, read_name.c_str())
      .def(py::init<Write<Scalar>>())
      .def("size", &Read<Scalar>::size);
  py::class_<Wrapper, Read<Scalar>>(module, py_wrapper.c_str())
      .def(py::init<Write<Scalar>>())
      .def(py::init<LO, Scalar, std::string const&>(), py::arg("size"),
          py::arg("value"), py::arg("name") = "");
  py::class_<HostRead<Scalar>>(
      module, hostread_name.c_str(), py::buffer_protocol())
      .def(py::init<Read<Scalar>>())
      .def_buffer([](HostRead<Scalar>& a) -> py::buffer_info {
        return py::buffer_info(
            /* note: although this breaks const-correctness, I think
               the benefits may outweigh the drawbacks */
            const_cast<Scalar*>(a.data()), sizeof(Scalar),
            py::format_descriptor<Scalar>::format(), 1, {a.size()},
            {sizeof(Scalar)});
      })
      .def("get", &HostRead<Scalar>::get);
  py::class_<HostWrite<Scalar>>(
      module, hostwrite_name.c_str(), py::buffer_protocol())
      .def(py::init<LO, std::string const&>(), py::arg("size"),
          py::arg("name") = "")
      .def_buffer([](HostWrite<Scalar>& a) -> py::buffer_info {
        return py::buffer_info(a.data(), sizeof(Scalar),
            py::format_descriptor<Scalar>::format(), 1, {a.size()},
            {sizeof(Scalar)});
      })
      .def("set", &HostWrite<Scalar>::set)
      .def("get", &HostWrite<Scalar>::get)
      .def("write", &HostWrite<Scalar>::write);
  Write<Scalar> (*deep_copy_type)(Read<Scalar> a, std::string const&) =
      &deep_copy;
  module.def(deepcopy_name.c_str(), deep_copy_type, py::arg("a"),
      py::arg("name") = "");
}

void pybind11_array(py::module& module) {
  pybind11_array_type<I8, Bytes>(module, "int8", "Bytes");
  pybind11_array_type<I32, LOs>(module, "int32", "LOs");
  pybind11_array_type<I64, GOs>(module, "int64", "GOs");
  pybind11_array_type<Real, Reals>(module, "float64", "Reals");
}

}  // namespace Omega_h
