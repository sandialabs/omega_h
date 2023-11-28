#include <Omega_h_file.hpp>


int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());

    int dim = mesh.dim();
    for(int d=0; d < mesh.dim(); d++)
        printf("(Dim, Num Mesh Entities): (%d, %d)\n", dim, mesh.nents(d));
}