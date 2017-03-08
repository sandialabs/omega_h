#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "Omega_h_array_ops.hpp"
#include "host_few.hpp"

namespace Omega_h {

class Mesh;

constexpr Int nhistogram_buckets = 10;

template <Int n>
struct Histogram {
  Real min;
  Real max;
  HostFew<GO, n> counts;
};

template <Int n>
Histogram<n> get_histogram(
    Mesh* mesh, Int dim, Reals values, Real min_value, Real max_value);
template <Int n>
void print_histogram(
    Mesh* mesh, Histogram<n> histogram, std::string const& name);

void print_goal_stats(Mesh* mesh, char const* name, Int ent_dim, Reals values,
    MinMax<Real> desired, MinMax<Real> actual);

#define INST_DECL(n)                                                           \
  extern template Histogram<n> get_histogram<n>(                               \
      Mesh * mesh, Int dim, Reals values, Real min_value, Real max_value);     \
  extern template void print_histogram<n>(                                     \
      Mesh * mesh, Histogram<n> histogram, std::string const& name);
INST_DECL(10)
#undef INST_DECL
}  // namespace Omega_h

#endif
