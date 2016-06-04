void modify_conn(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    LOs prod_verts2verts,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs old_lows2new_lows);

void modify_owners(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs old_ents2new_ents);

LOs collect_same(Mesh& mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds);

LOs get_keys2reps(Mesh& mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2kds,
    LOs keys2nprods);

LOs get_rep_counts(
    Mesh& mesh,
    Int ent_dim,
    LOs keys2reps,
    LOs keys2nprods,
    LOs same_ents2ents);

template <typename T>
void find_new_offsets(
    Int ent_dim,
    Read<T> old_ents2new_offsets,
    LOs same_ents2old_ents,
    LOs keys2reps,
    LOs keys2prods,
    Read<T>& same_ents2new_offsets,
    Read<T>& prods2new_offsets);

void modify_globals(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2ents,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents);
