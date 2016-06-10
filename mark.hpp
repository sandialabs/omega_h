Read<I8> mark_exposed_sides(Mesh& mesh);
Read<I8> mark_down(Mesh& mesh, Int high_dim, Int low_dim,
    Read<I8> marked_highs);
Read<I8> mark_up(Mesh& mesh, Int low_dim, Int high_dim,
    Read<I8> low_marked);
Read<I8> mark_by_class_dim(Mesh& mesh, Int ent_dim, Int class_dim);
Read<I8> mark_by_owner(Mesh& mesh, Int ent_dim, I32 rank);
