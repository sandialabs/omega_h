/* if the mesh satisfies quality and
   length requirements, print a short
   message and return true.
   otherwise, print a more detailed
   statistical report and return false */
bool adapt_check(Mesh* mesh,
    Real qual_floor,
    Real qual_ceil,
    Real len_floor,
    Real len_ceil);

/* returns true if the mesh topology was modified. */
bool adapt(Mesh* mesh,
    Real qual_ceil,
    Real qual_floor,
    Real len_floor,
    Real len_ceil,
    Int nlayers);
