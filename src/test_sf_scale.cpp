#include <iostream>

#include "omega_h.hpp"

#include "size.hpp"

int main(int argc, char** argv) {
  auto lib = osh::Library(&argc, &argv);
  osh::Mesh mesh;
  osh::build_box(&mesh, lib, 1, 1, 1, 8, 8, 8);
  auto target_nelems = mesh.nelems();
  {
    auto size = osh::find_identity_size(&mesh);
    auto elems_per_elem = expected_elems_per_elem_iso(&mesh, size);
    auto elems = repro_sum_owned(&mesh, mesh.dim(), elems_per_elem);
    auto size_scal = target_nelems / elems;
    std::cout << "ratio for iso " << size_scal << '\n';
  }
  {
    auto metric = osh::find_identity_metric(&mesh);
    auto elems_per_elem = expected_elems_per_elem_metric(&mesh, metric);
    auto elems = repro_sum_owned(&mesh, mesh.dim(), elems_per_elem);
    auto size_scal = target_nelems / elems;
    std::cout << "ratio for aniso " << size_scal << '\n';
  }
}
