namespace inertia {

Read<I8> mark_bisection(CommPtr comm,
    Reals coords, Reals masses, Real tolerance);
Read<I8> mark_bisection(CommPtr comm,
    Reals coords, Reals masses, Real tolerance,
    Vector<3> axis);

}
