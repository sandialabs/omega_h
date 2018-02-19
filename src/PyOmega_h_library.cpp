#include <Omega_h_library.hpp>
#include <PyOmega_h.hpp>
#include <iostream>

namespace Omega_h {

/* The lifetime of the Library object is quite important (it must contain the
   lifetime of pretty much all other Omega_h objects), and
   I'm unsure about the order in which destructors will be called.
   So, for now, I'll take the approach that the Python interface will have
   the Library as a hidden global variable.
   This is consistent with how mpi4py seems to work.
   I tried using the Python atexit mechanism, but that seems to execute prior
   to final garbage collection.
 */

std::unique_ptr<Library> pybind11_global_library;

void pybind11_library(py::module& module) {
  pybind11_global_library =
      decltype(pybind11_global_library)(new Omega_h::Library());
  module.def("world", []() { return pybind11_global_library->world(); },
      "Returns the world communicator");
  module.def("self", []() { return pybind11_global_library->self(); },
      "Returns the self communicator");
}

}  // namespace Omega_h
