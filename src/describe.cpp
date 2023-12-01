#include <Omega_h_file.hpp>


int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());
    std::ostringstream oss;

    #ifdef OMEGA_H_USE_MPI
        int comm_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
        oss << "\nComm Rank: " << comm_rank << "\n";
    #endif

    oss << "\nMesh Entity Count by Dimension: (Dim, Entities Count)\n";
    for(int dim=0; dim < mesh.dim(); dim++)
        oss << "(" << dim << ", " << mesh.nents(dim) << ")\n";

    oss << "\nImbalance by Dimension: (Dim, Imbalance)\n";
    for(int dim=0; dim < mesh.dim(); dim++)
        oss << "(" << dim << ", " << mesh.imbalance(dim) << ")\n";

    // oss << "\nShapes:\n";
    // oss << "Num Pyrams: " << mesh.npyrams() << "\n";
    // oss << "Num Wedges: " << mesh.nwedges() << "\n";
    // oss << "Num Hexs: " << mesh.nhexs() << "\n";
    // oss << "Num Tets: " << mesh.ntets() << "\n";
    // oss << "Num Quads: " << mesh.nquads() << "\n";
    // oss << "Num Tris: " << mesh.ntris() << "\n";

    oss << "\nTags by Dimension: (Dim, Tag, Size per Entity)\n";
    for (int dim=0; dim < mesh.dim(); dim++)
    for (int tag=0; tag < mesh.ntags(dim); tag++) {
        auto tagbase = mesh.get_tag(dim, tag);
        int size;
        if (tagbase->type() == OMEGA_H_I8)
            size = mesh.get_tag<Omega_h::I8>(dim, tagbase->name())->array().size();
        if (tagbase->type() == OMEGA_H_I32)
            size = mesh.get_tag<Omega_h::I32>(dim, tagbase->name())->array().size();
        if (tagbase->type() == OMEGA_H_I64)
            size = mesh.get_tag<Omega_h::I64>(dim, tagbase->name())->array().size();
        if (tagbase->type() == OMEGA_H_F64)
            size = mesh.get_tag<Omega_h::Real>(dim, tagbase->name())->array().size();

        size /= mesh.nents(dim);
        oss << "(" << dim << ", " << tagbase->name().c_str() << ", " << size << ")\n";
    }

    oss << "\n--------------------------------------------------------\n";
    std::cout << oss.str();
}