namespace surf {
Reals get_side_normals(Mesh& mesh, LOs surf_side2side);
Reals get_hinge_angles(Mesh& mesh,
    Reals surf_side_normals,
    LOs surf_hinge2hinge,
    LOs side2surf_side);
}
