#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include <vector>

#include "Omega_h_array_ops.hpp"
#include "host_few.hpp"

namespace Omega_h {

class Mesh;

struct Histogram {
  Real min;
  Real max;
  std::vector<GO> bins;
};

Histogram get_histogram(
    Mesh* mesh, Int dim, Int nbins, Real min_value, Real max_value, Reals values);

void print_histogram(
    Mesh* mesh, Histogram const& histogram, std::string const& name);

void print_goal_stats(Mesh* mesh, char const* name, Int ent_dim, Reals values,
    MinMax<Real> desired, MinMax<Real> actual);

}  // namespace Omega_h

#endif
