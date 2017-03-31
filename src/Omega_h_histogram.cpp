#include "Omega_h_histogram.hpp"

#include <iomanip>
#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

Histogram get_histogram(
    Mesh* mesh, Int dim, Int nbins, Real min_value, Real max_value, Reals values) {
  auto owned_values = mesh->owned_array(dim, values, 1);
  auto interval = (max_value - min_value) / Real(nbins);
  Histogram histogram;
  histogram.min = min_value;
  histogram.max = max_value;
  histogram.bins.resize(std::size_t(nbins));
  for (Int i = 0; i < nbins; ++i) {
    auto floor = interval * i + min_value;
    auto ceil = interval * (i + 1) + min_value;
    auto marked = land_each(
        each_geq_to(owned_values, floor), each_lt(owned_values, ceil));
    histogram.bins[std::size_t(i)] = get_sum(mesh->comm(), marked);
  }
  return histogram;
}

void print_histogram(
    Mesh* mesh, Histogram const& histogram, std::string const& name) {
  if (mesh->comm()->rank()) return;
  std::ios saved_state(0);
  saved_state.copyfmt(std::cout);
  std::cout << std::fixed << std::setprecision(2);
  std::cout << name << " histogram:\n";
  auto nbins = Int(histogram.bins.size());
  auto interval = (histogram.max - histogram.min) / Real(nbins);
  for (Int i = 0; i < nbins; ++i) {
    auto floor = interval * i + histogram.min;
    auto ceil = interval * (i + 1) + histogram.min;
    std::cout << floor << '-' << ceil << ": "
      << histogram.bins[std::size_t(i)] << '\n';
  }
  std::cout.copyfmt(saved_state);
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

}  // namespace Omega_h
