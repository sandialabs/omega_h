#include <Omega_h_element.hpp>
#include <Omega_h_file.hpp>


int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    auto comm = lib.world();
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());

    auto verbose = false;
    if (argc == 3) verbose = (std::string(argv[2]) == "on");

    const int rank = comm->rank();

    std::array<Omega_h::GO, 4> counts;
    for(int dim=0; dim < mesh.dim(); dim++)
       counts[dim] = mesh.nglobal_ents(dim);

    std::array<double, 4> imb;
    for(int dim=0; dim < mesh.dim(); dim++)
       imb[dim] = mesh.imbalance(dim);

    std::ostringstream oss;
    // always print two places to the right of the decimal
    // for floating point types (i.e., imbalance)
    oss.precision(2);
    oss << std::fixed;

    if(!rank) {
        oss << "\nMesh Entity Type: " << Omega_h::topological_singular_name(mesh.family(), mesh.dim()-1) << "\n";

        oss << "\nGlobal Mesh Entity Count and Imbalance (max/avg): (Dim, Entity Count, Imbalance)\n";
        for(int dim=0; dim < mesh.dim(); dim++)
            oss << "(" << dim << ", " << counts[dim] << ", " << imb[dim] << ")\n";

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

        std::cout << oss.str();
    }

    if(verbose) {
        comm->barrier(); // write the per-part data at the end
        if(!rank) {
          std::cout << "\nPer Rank Mesh Entity Count: (Rank: Entity Count by Dim <0,1,2,3>)\n";
        }
        comm->barrier(); // write the per-part data at the end
        oss.str(""); // clear the stream

        std::array<Omega_h::LO, 4> counts = {0,0,0,0};
        for(int dim=0; dim < mesh.dim(); dim++)
            counts[dim] = mesh.nents(dim);
        oss << "(" << rank << ": " << counts[0] << ", "
                                        << counts[1] << ", "
                                        << counts[2] << ", "
                                        << counts[3] << ")\n";
        std::cout << oss.str();
    }
    return 0;
}
