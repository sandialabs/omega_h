#include <Omega_h_file.hpp>


int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());

    int dim = mesh.dim();
    for(int d=0; d < mesh.dim(); d++)
        printf("(Dim, Num Mesh Entities): (%d, %d)\n", d, mesh.nents(d));

    printf("nelems: %d\n", mesh.nelems());
    printf("nregions: %d\n", mesh.nregions());
    printf("nfaces: %d\n", mesh.nfaces());
    printf("nedges: %d\n", mesh.nedges());
    printf("nverts: %d\n", mesh.nverts());

    printf("npyrams: %d\n", mesh.npyrams());
    printf("nwedges: %d\n", mesh.nwedges());
    printf("nhexs: %d\n", mesh.nhexs());
    printf("ntets: %d\n", mesh.ntets());
    printf("nquads: %d\n", mesh.nquads());
    printf("ntris: %d\n", mesh.ntris());
}