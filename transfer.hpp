void transfer_linear_interp(Mesh& old_mesh, Mesh& new_mesh,
    LOs key_verts2verts,
    LOs keys2midverts,
    LOs same_verts2old_verts,
    LOs same_verts2new_verts);

void transfer_metric(Mesh& old_mesh, Mesh& new_mesh,
    LOs key_verts2verts,
    LOs keys2midverts,
    LOs same_verts2old_verts,
    LOs same_verts2new_verts);

void transfer_inherit_refine(Mesh& old_mesh, Mesh& new_mesh,
    LOs keys2edges,
    Int prod_dim,
    LOs keys2prods,
    LOs prods2new_ents,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents);

void transfer_length(Mesh& old_mesh, Mesh& new_mesh,
    LOs same_ents2old_ents,
    LOs same_ents2new_ents,
    LOs prods2new_ents);
