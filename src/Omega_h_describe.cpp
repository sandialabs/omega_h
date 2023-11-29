#include <Omega_h_file.hpp>


int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());

    printf("\nMesh Entity Count by Dimension: (Dim, Num Mesh Entities)\n");
    for(int dim=0; dim < mesh.dim(); dim++)
        printf("(%d, %d)\n", dim, mesh.nents(dim));

    printf("\nImbalance by Dimension: (Dim, Imbalance)\n");
    for(int dim=0; dim < mesh.dim(); dim++)
        printf("(%d, %d)\n", dim, mesh.imbalance(dim));

    printf("\nShapes:\n");
    printf("Num Pyrams: %d\n", mesh.npyrams());
    printf("Num Wedges: %d\n", mesh.nwedges());
    printf("Num Hexs: %d\n", mesh.nhexs());
    printf("Num Tets: %d\n", mesh.ntets());
    printf("Num Quads: %d\n", mesh.nquads());
    printf("Num Tris: %d\n", mesh.ntris());

    printf("\nTags: (Dim, Tag, Size)\n");
    for (int dim=0; dim < mesh.dim(); dim++)
    for (int tag=0; tag < mesh.ntags(dim); tag++) {
        auto tagbase = mesh.get_tag(dim, tag);
        printf("(%d, %s, ", dim, tagbase->name().c_str());

        if (tagbase->type() == OMEGA_H_I8)
            printf("%d)\n", mesh.get_tag<Omega_h::I8>(dim, tagbase->name())->array().size());
        if (tagbase->type() == OMEGA_H_I32)
            printf("%d)\n", mesh.get_tag<Omega_h::I32>(dim, tagbase->name())->array().size());
        if (tagbase->type() == OMEGA_H_I64)
            printf("%d)\n", mesh.get_tag<Omega_h::I64>(dim, tagbase->name())->array().size());
        if (tagbase->type() == OMEGA_H_F64)
            printf("%d)\n", mesh.get_tag<Omega_h::Real>(dim, tagbase->name())->array().size());
    }
}