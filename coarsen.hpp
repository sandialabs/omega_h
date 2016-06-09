Few<Read<I8>, 4> mark_dead_ents(Mesh& mesh,
    LOs rails2edges,
    Read<I8> rail_col_verts);

Adj find_coarsen_domains(Mesh& mesh,
    LOs keys2verts,
    Int ent_dim,
    Read<I8> ents_are_dead);

Reals coarsen_qualities(Mesh& mesh,
    LOs cands2edges,
    Read<I8> cand_codes);

void choose_vertex_collapses(Mesh& mesh,
    LOs cands2edges,
    Read<I8> cand_edge_codes,
    Reals cand_edge_quals,
    Read<I8>& verts_are_cands,
    Reals& vert_quals);

LOs coarsen_topology(Mesh& mesh,
    LOs keys2verts_onto,
    Int dom_dim,
    Adj keys2doms);
