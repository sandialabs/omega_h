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

void modify_globals(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim,
    Int key_dim,
    LOs keys2ents,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents);
