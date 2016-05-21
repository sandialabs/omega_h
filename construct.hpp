void build_from_elems2verts(Mesh& mesh, Int edim, LOs ev2v, LO nverts);
void build_from_elems_and_coords(Mesh& mesh, Int edim, LOs ev2v, Reals coords);
void build_box(Mesh& mesh,
    Real x, Real y, Real z,
    LO nx, LO ny, LO nz);
