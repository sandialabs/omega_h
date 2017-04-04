#include "Omega_h_histogram.hpp"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Omega_h_array_ops.hpp"
#include "Omega_h_simplex.hpp"

namespace Omega_h {

Histogram get_histogram(Mesh* mesh, Int dim, Int nbins, Real min_value,
    Real max_value, Reals values) {
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

void print_histogram(Histogram const& histogram, std::string const& name) {
  std::ios saved_state(0);
  saved_state.copyfmt(std::cout);
  std::cout << std::fixed << std::setprecision(2);
  std::cout << name << " histogram:\n";
  auto nbins = Int(histogram.bins.size());
  auto interval = (histogram.max - histogram.min) / Real(nbins);
  for (Int i = 0; i < nbins; ++i) {
    auto floor = interval * i + histogram.min;
    auto ceil = interval * (i + 1) + histogram.min;
    std::cout << floor << '-' << ceil << ": " << histogram.bins[std::size_t(i)]
              << '\n';
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

void render_histogram_matplotlib(Histogram const& histogram,
    std::string const& filepath, std::string const& ents_name,
    std::string const& var_name) {
  std::ofstream script("Omega_h_histogram.py");
  script << "#!/usr/bin/python\n";
  script << "import matplotlib\n";
  script << "# Force matplotlib to not use any Xwindows backend.\n";
  script << "matplotlib.use('Agg')\n";
  script << "import matplotlib.pyplot as plt\n";
  auto nbins = Int(histogram.bins.size());
  auto interval = (histogram.max - histogram.min) / Real(nbins);
  script << "a = [";
  GO binmax = 0;
  for (Int i = 0; i < nbins; ++i) {
    auto mid = (interval * i) + histogram.min + (interval / 2.0);
    script << mid;
    if (i + 1 < nbins) script << ", ";
    binmax = max2(binmax, histogram.bins[std::size_t(i)]);
  }
  script << "]\n";
  script << "b = [";
  for (Int i = 0; i < nbins; ++i) {
    script << histogram.bins[std::size_t(i)];
    if (i + 1 < nbins) script << ", ";
  }
  script << "]\n";
  script << "ax = plt.gca()\n";
  script << "ax.set_autoscale_on(False)\n";
  script << "plt.hist(a, " << nbins << ", weights=b";
  script << ", range=(" << histogram.min << ", " << histogram.max << ")";
  script << ")\n";
  script << "plt.title('Histogram of " << ents_name << " by " << var_name
         << "')\n";
  script << "plt.axis([" << histogram.min << ", " << histogram.max << ", ";
  script << "0, " << binmax << "])\n";
  script << "plt.savefig('" << filepath << "', bbox_inches='tight')\n";
  script.close();
  int ret = ::system("python Omega_h_histogram.py");
  OMEGA_H_CHECK(ret == 0);
}

}  // namespace Omega_h
