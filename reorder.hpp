Graph find_entities_of_first_vertices(
    Mesh& mesh, Int ent_dim);
LOs ent_order_from_vert_order(Mesh& mesh,
    Int ent_dim, LOs new_vert2old_vert);
void reorder_tags(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim, LOs new_ents2old_ents);
void reorder_down(Mesh& old_mesh, Mesh& new_mesh,
    Int ent_dim, LOs new_ents2old_ents,
    LOs old_lows2new_lows);
void reorder_mesh(Mesh& old_mesh, Mesh& new_mesh,
    LOs new_verts2old_verts);
void reorder_mesh(Mesh& mesh,
    LOs new_verts2old_verts);
void reorder_by_hilbert(Mesh& mesh);
