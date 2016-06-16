void swap3d_qualities(Mesh* mesh, LOs cands2edges,
    Reals* cand_quals, Read<I8>* cand_configs);

std::array<LOs, 4> swap3d_keys_to_prods(Mesh* mesh, LOs keys2edges);

std::array<LOs, 4> swap3d_topology(Mesh* mesh,
    LOs keys2edges,
    Read<I8> edge_configs,
    std::array<LOs, 4> keys2prods);

bool run_swap3d(Mesh* mesh, Real qual_ceil, Int nlayers);
