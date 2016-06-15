Reals swap2d_qualities(Mesh* mesh, LOs cands2edges);
void filter_swap2d_improve(Mesh* mesh,
    LOs* cands2edges, Reals* cand_quals);
void swap2d_topology(Mesh* mesh, LOs keys2edges,
    LOs* prod_tv2v, LOs* prod_ev2v);
