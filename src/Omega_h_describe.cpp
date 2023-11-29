#include <Omega_h_file.hpp>


int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());

    printf("\nMesh Properties:\n");
    for(int dim=0; dim < mesh.dim(); dim++)
        printf("(Dim, Num Mesh Entities): (%d, %d)\n", dim, mesh.nents(dim));

    printf("\n");
    for(int dim=0; dim < mesh.dim(); dim++)
        printf("(Dim, Imbalance): (%d, %d)\n", dim, mesh.imbalance(dim));

    printf("\nGeometry:\n");
    printf("Num Elems: %d\n", mesh.nelems());
    printf("Num Regions: %d\n", mesh.nregions());
    printf("Num Faces: %d\n", mesh.nfaces());
    printf("Num Edges: %d\n", mesh.nedges());
    printf("Num Verts: %d\n", mesh.nverts());

    printf("\nShapes:\n");
    printf("Num Pyrams: %d\n", mesh.npyrams());
    printf("Num Wedges: %d\n", mesh.nwedges());
    printf("Num Hexs: %d\n", mesh.nhexs());
    printf("Num Tets: %d\n", mesh.ntets());
    printf("Num Quads: %d\n", mesh.nquads());
    printf("Num Tris: %d\n", mesh.ntris());

    printf("\nTags:\n");
    for (int dim=0; dim < mesh.dim(); dim++)
    for (int tag=0; tag < mesh.ntags(dim); tag++) {
        auto tagbase = mesh.get_tag(dim, tag);
        printf("(Dim, Tag): (%d, %s)\n", dim, tagbase->name().c_str());
    }
}