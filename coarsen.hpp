std::array<Read<I8>, 4> mark_dead_ents(Mesh* mesh, LOs rails2edges,
                                       Read<I8> rail_col_dirs);

Adj find_coarsen_domains(Mesh* mesh, LOs keys2verts, Int ent_dim,
                         Read<I8> ents_are_dead);

Reals coarsen_qualities(Mesh* mesh, LOs cands2edges, Read<I8> cand_codes);

Read<I8> filter_coarsen_min_qual(Read<I8> cand_codes, Reals cand_quals,
                                 Real min_qual);

Read<I8> filter_coarsen_improve(Mesh* mesh, LOs cands2edges,
                                Read<I8> cand_codes, Reals cand_quals);

void choose_vertex_collapses(Mesh* mesh, LOs cands2edges,
                             Read<I8> cand_edge_codes, Reals cand_edge_quals,
                             Read<I8>& verts_are_cands, Reals& vert_quals);

void find_rails(Mesh* mesh, LOs keys2verts, Reals vert_quals,
                Read<I8> edge_cand_codes, Reals edge_cand_quals,
                LOs& rails2edges, Read<I8>& rail_col_dirs);

LOs get_verts_onto(Mesh* mesh, LOs rails2edges, Read<I8> rail_col_dirs);

LOs coarsen_topology(Mesh* mesh, LOs keys2verts_onto, Int dom_dim,
                     Adj keys2doms, LOs old_verts2new_verts);

bool coarsen(Mesh* mesh, Real min_qual, bool improve, bool verbose);

bool coarsen_verts(Mesh* mesh, Read<I8> vert_marks, Real min_qual, bool improve,
                   bool verbose);

bool coarsen_ents(Mesh* mesh, Int ent_dim, Read<I8> marks, Real min_qual,
                  bool improve, bool verbose);

bool coarsen_by_size(Mesh* mesh, Real min_len, Real min_qual, bool verbose);

bool coarsen_slivers(Mesh* mesh, Real qual_ceil, Int nlayers, bool verbose);
