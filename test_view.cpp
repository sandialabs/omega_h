#include "omega_h.hpp"

int main(int argc, char** argv) {
//this caused a Kokkos::initialize segfault...
//auto library = osh::Library(nullptr, nullptr);
  auto library = osh::Library(&argc, &argv);
  osh::Mesh mesh;
  osh::binary::read(argv[1], library.world(), &mesh);
  mesh.balance();
  mesh.reorder();
  std::vector<osh::Int> class_dims;
  std::vector<osh::I32> class_ids;
  class_dims.push_back(2);
  class_ids.push_back(85);
  osh::Int dim = 0;
  auto nents = mesh.nents(dim);
  osh::Read<osh::I8> marked_read = osh::mark_class_closures(&mesh, dim,
      class_dims, class_ids);
  assert(marked_read.size() == nents);
  std::cerr << "marked mesh entities\n";
  Kokkos::View<osh::I8 const*> marked = marked_read.view();
  OSH_CHECK(marked.size() == nents);
  std::cerr << "got marked = marked_read.view()\n";
  typename Kokkos::View<osh::I8 const*>::HostMirror mirror = Kokkos::create_mirror_view(marked);
  std::cerr << "ran mirror = Kokkos::create_mirror_view(marked)\n";
  Kokkos::deep_copy(mirror, marked);
  std::cerr << "ran Kokkos::deep_copy(mirror, marked)\n";
  for (decltype(nents) i = 0; i < nents; ++i)
    if (mirror(i))
      std::cerr << i << " in set\n";
}
