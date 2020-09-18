#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <cstdlib>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  OMEGA_H_CHECK(argc == 3 || argc == 4);
  Omega_h::Mesh mesh(&lib);
  printf("osh2vtk ok1\n");
  Omega_h::binary::read(argv[1], lib.world(), &mesh);
  printf("osh2vtk ok2\n");
  auto dim = mesh.dim();
  if (argc == 4) dim = atoi(argv[2]);
  printf("osh2vtk ok3\n");
  Omega_h::vtk::write_parallel(argv[argc - 1], &mesh, dim);
  printf("osh2vtk ok4\n");
}
