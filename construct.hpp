void add_ents2verts(Mesh* mesh, Int edim, LOs ev2v,
    Read<GO> vert_globals);
void resolve_derived_copies(
    CommPtr comm,
    Read<GO> verts2globs,
    Int deg,
    LOs* p_ent_verts2verts,
    Remotes* p_ents2owners);
