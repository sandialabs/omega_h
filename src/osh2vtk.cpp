#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include <cstdlib>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  if( argc != 3 && argc != 4 ) {
    fprintf(stderr, "Usage: %s inputMesh.osh [dimension=0|1|2|3] outputMesh.vtk\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 3 || argc == 4);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], lib.world(), &mesh);
  auto dim = mesh.dim();
  if (argc == 4) dim = atoi(argv[2]);
  OMEGA_H_CHECK(dim>=0 && dim<=mesh.dim());
  Omega_h::vtk::write_parallel(argv[argc - 1], &mesh, dim);
}
