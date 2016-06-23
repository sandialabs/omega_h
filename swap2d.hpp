Reals swap2d_qualities(Mesh* mesh, LOs cands2edges);

void swap2d_topology(Mesh* mesh, LOs keys2edges,
    std::array<LOs, 3>* keys2prods,
    std::array<LOs, 3>* prod_verts2verts);

bool swap2d(Mesh* mesh, Real qual_ceil, Int nlayers, bool verbose);
