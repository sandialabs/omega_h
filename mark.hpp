Read<I8> mark_exposed_sides(Mesh& mesh);
Read<I8> mark_down(Mesh& mesh, Int high_dim, Int low_dim,
    Read<I8> marked_highs);

template <typename T>
Read<I8> mark_equal(Read<T> a, T val);

Read<I8> mark_by_class_dim(Mesh& mesh, Int ent_dim, Int class_dim);
Read<I8> mark_by_owner(Mesh& mesh, Int ent_dim, I32 rank);
