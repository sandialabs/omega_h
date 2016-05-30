namespace inertia {

Read<I8> mark_bisection(CommPtr comm,
    Reals coords, Reals masses, Real tolerance,
    Vector<3>& axis);
Read<I8> mark_bisection_given_axis(CommPtr comm,
    Reals coords, Reals masses, Real tolerance,
    Vector<3> axis);

}
