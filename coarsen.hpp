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
