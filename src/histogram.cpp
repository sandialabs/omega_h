#include "histogram.hpp"
#include <iomanip>
#include <iostream>
#include "Omega_h_array_ops.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

template <Int n>
Histogram<n> get_histogram(
    Mesh* mesh, Int dim, Reals values, Real min_value, Real max_value) {
  auto owned_values = mesh->owned_array(dim, values, 1);
  auto interval = Real(max_value - min_value) / n;
  Histogram<n> histogram;
  histogram.min = min_value;
  histogram.max = max_value;
  for (Int i = 0; i < n; ++i) {
    auto floor = interval * i + min_value;
    auto ceil = interval * (i + 1) + min_value;
    auto marked = land_each(
        each_geq_to(owned_values, floor), each_lt(owned_values, ceil));
    histogram.counts[i] = get_sum(mesh->comm(), marked);
  }
  return histogram;
}

template <Int n>
void print_histogram(
    Mesh* mesh, Histogram<n> histogram, std::string const& name) {
  if (mesh->comm()->rank()) return;
  auto precision_before = std::cout.precision();
  std::ios::fmtflags stream_state(std::cout.flags());
  std::cout << std::fixed << std::setprecision(2);
  std::cout << name << " histogram:\n";
  auto interval = (histogram.max - histogram.min) / n;
  for (Int i = 0; i < n; ++i) {
    auto floor = interval * i + histogram.min;
    auto ceil = interval * (i + 1) + histogram.min;
    std::cout << floor << '-' << ceil << ": " << histogram.counts[i] << '\n';
  }
  std::cout.flags(stream_state);
  std::cout.precision(precision_before);
}

void print_goal_stats(Mesh* mesh, char const* name, Int ent_dim, Reals values,
    MinMax<Real> desired, MinMax<Real> actual) {
  auto low_marks = each_lt(values, desired.min);
  auto high_marks = each_gt(values, desired.max);
  auto nlow = count_owned_marks(mesh, ent_dim, low_marks);
  auto nhigh = count_owned_marks(mesh, ent_dim, high_marks);
  auto ntotal = mesh->nglobal_ents(ent_dim);
  auto nmid = ntotal - nlow - nhigh;
  if (mesh->comm()->rank() == 0) {
    auto precision_before = std::cout.precision();
    std::ios::fmtflags stream_state(std::cout.flags());
    std::cout << std::fixed << std::setprecision(2);
    std::cout << ntotal << " " << plural_names[ent_dim];
    std::cout << ", " << name << " [" << actual.min << "," << actual.max << "]";
    if (nlow) {
      std::cout << ", " << nlow << " <" << desired.min;
    }
    if (nmid) {
      std::cout << ", " << nmid << " in [" << desired.min << "," << desired.max
                << "]";
    }
    if (nhigh) {
      std::cout << ", " << nhigh << " >" << desired.max;
    }
    std::cout << '\n';
    std::cout.flags(stream_state);
    std::cout.precision(precision_before);
  }
}

#define INST(n)                                                                \
  template Histogram<n> get_histogram<n>(                                      \
      Mesh * mesh, Int dim, Reals values, Real min_value, Real max_value);     \
  template void print_histogram<n>(                                            \
      Mesh * mesh, Histogram<n> histogram, std::string const& name);
INST(10)
#undef INST
}  // namespace Omega_h
