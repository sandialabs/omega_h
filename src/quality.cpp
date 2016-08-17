#include "quality.hpp"

#include <iomanip>
#include <iostream>

#include "array.hpp"
#include "loop.hpp"
#include "map.hpp"

namespace Omega_h {

template <typename ElementQualities, Int dim>
Reals measure_qualities_tmpl(Mesh* mesh, LOs a2e) {
  ElementQualities measurer(mesh);
  auto ev2v = mesh->ask_verts_of(mesh->dim());
  auto na = a2e.size();
  Write<Real> qualities(na);
  auto f = LAMBDA(LO a) {
    auto e = a2e[a];
    auto v = gather_verts<dim + 1>(ev2v, e);
    qualities[a] = measurer.measure(v);
  };
  parallel_for(na, f);
  return qualities;
}

Reals measure_qualities(Mesh* mesh, LOs a2e) {
  if (mesh->dim() == 3) {
    if (mesh->has_tag(VERT, "metric")) {
      return measure_qualities_tmpl<MetricElementQualities, 3>(mesh, a2e);
    } else {
      return measure_qualities_tmpl<RealElementQualities, 3>(mesh, a2e);
    }
  } else {
    CHECK(mesh->dim() == 2);
    if (mesh->has_tag(VERT, "metric")) {
      return measure_qualities_tmpl<MetricElementQualities, 2>(mesh, a2e);
    } else {
      return measure_qualities_tmpl<RealElementQualities, 2>(mesh, a2e);
    }
  }
}

Reals measure_qualities(Mesh* mesh) {
  return measure_qualities(mesh, LOs(mesh->nelems(), 0, 1));
}

QualityHistogram get_quality_histogram(Mesh* mesh) {
  Reals owned_qualities;
  if (mesh->could_be_shared(mesh->dim())) {
    auto overlap_qualities = mesh->ask_qualities();
    auto owned = mesh->owned(mesh->dim());
    auto owned2overlap = collect_marked(owned);
    owned_qualities = unmap(owned2overlap, overlap_qualities, 1);
  } else {
    owned_qualities = mesh->ask_qualities();
  }
  QualityHistogram histogram;
  auto interval = Real(1.0) / nquality_histogram_buckets;
  for (Int i = 0; i < nquality_histogram_buckets; ++i) {
    auto floor = interval * i;
    auto ceil = interval * (i + 1);
    auto marked = land_each(
        each_geq_to(owned_qualities, floor), each_lt(owned_qualities, ceil));
    auto nlocal_marked = sum(marked);
    auto nglobal_marked = mesh->comm()->allreduce(nlocal_marked, OMEGA_H_SUM);
    histogram[i] = nglobal_marked;
  }
  return histogram;
}

void print_quality_histogram(QualityHistogram histogram) {
  auto precision_before = std::cout.precision();
  std::ios::fmtflags stream_state(std::cout.flags());
  std::cout << std::fixed << std::setprecision(1);
  std::cout << "quality histogram:\n";
  auto interval = Real(1.0) / nquality_histogram_buckets;
  for (Int i = 0; i < nquality_histogram_buckets; ++i) {
    auto floor = interval * i;
    auto ceil = interval * (i + 1);
    std::cout << floor << '-' << ceil << ": " << histogram[i] << '\n';
  }
  std::cout.flags(stream_state);
  std::cout.precision(precision_before);
}

}  // end namespace Omega_h
