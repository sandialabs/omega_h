#ifndef INERTIA_HPP
#define INERTIA_HPP

#include "algebra.hpp"

namespace osh {

namespace inertia {

struct Rib {
  std::vector<Vector<3>> axes;
};

Read<I8> mark_bisection(
    CommPtr comm, Reals coords, Reals masses, Real tolerance, Vector<3>& axis);
Read<I8> mark_bisection_given_axis(
    CommPtr comm, Reals coords, Reals masses, Real tolerance, Vector<3> axis);
Rib recursively_bisect(CommPtr comm, Reals& coords, Reals& masses,
    Remotes& owners, Real tolerance, Rib hints);
}

}  // end namespace osh

#endif
