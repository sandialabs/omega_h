#include <Omega_h_file.hpp>


int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());

    for(int dim=0; dim < mesh.dim(); dim++)
        printf("(Dim, Num Mesh Entities): (%d, %d)\n", dim, mesh.nents(dim));

    for(int dim=0; dim < mesh.dim(); dim++)
        printf("(Dim, Imbalance): (%d, %d)\n", dim, mesh.imbalance(dim));

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

    for(int dim=0; dim < mesh.dim(); dim++)
    for (int tag=0; tag < mesh.ntags(dim); tag++) {
        auto tagbase = mesh.get_tag(dim, tag);
        printf("(Dim, Tag): (%d, %s)\n", dim, tagbase->name().c_str());
    }
}