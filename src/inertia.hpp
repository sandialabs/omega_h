#ifndef INERTIA_HPP
#define INERTIA_HPP

#include <vector>

#include "Omega_h_remotes.hpp"
#include "Omega_h_vector.hpp"

namespace Omega_h {

namespace inertia {

struct Rib {
  std::vector<Vector<3>> axes;
};

Read<I8> mark_bisection(
    CommPtr comm, Reals coords, Reals masses, Real tolerance, Vector<3>& axis);
Read<I8> mark_bisection_given_axis(
    CommPtr comm, Reals coords, Reals masses, Real tolerance, Vector<3> axis);
void recursively_bisect(CommPtr comm, Real tolerance, Reals* p_coords,
    Reals* p_masses, Remotes* p_owners, Rib* p_hints);
}  // namespace inertia

}  // end namespace Omega_h

#endif
