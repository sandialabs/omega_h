#include <Omega_h_library.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_random.hpp>

using namespace Omega_h;

static OMEGA_H_DEVICE void contrib(Write<Real> const& buckets, Real value) {
  if (value < 0.0 || value >= 1.0) return;
  buckets[LO(std::floor(value * buckets.size()))]++;
}

static Real test_chi_squared(std::string const& name, HostRead<Real> buckets,
    HostRead<Real> cdfs, Real num_distribution_parameters, Real nsamples,
    Real cutoff) {
  OMEGA_H_CHECK(cdfs.size() == buckets.size() + 1);
  LO num_nonempty_buckets = 0;
  for (LO i = 0; i < buckets.size(); ++i) {
    if (buckets[i] != 0) ++num_nonempty_buckets;
  }
  Real chi_squared = 0;
  for (LO i = 0; i < buckets.size(); ++i) {
    auto expected = (cdfs[i + 1] - cdfs[i]) * nsamples;
    chi_squared += square(buckets[i] - expected) / expected;
  }
  auto ndofs = buckets.size() - num_distribution_parameters;
  auto p_value = 1.0 - cumulative_chi_squared_density(ndofs, chi_squared);
  if (p_value < cutoff) {
    Omega_h_fail(
        "Chi-squared test failed for distribution %s with p-value %f < %f\n",
        name.c_str(), p_value, cutoff);
  } else {
    std::cout << "Chi-squared test passed for distribution " << name
              << " with p-value " << p_value << " >= " << cutoff << '\n';
  }
  return p_value;
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  I64 seed = 1771;
  I64 key = 365;
  LO nbuckets = 100;
  int nsamples = 1000 * nbuckets;
  auto d_uniform_buckets = Write<Real>(nbuckets, 0);
  auto d_normal_buckets = Write<Real>(nbuckets, 0);
  auto d_weibull_1_buckets = Write<Real>(nbuckets, 0);
  auto d_weibull_2_buckets = Write<Real>(nbuckets, 0);
  auto d_weibull_3_buckets = Write<Real>(nbuckets, 0);
  auto d_weibull_4_buckets = Write<Real>(nbuckets, 0);
  auto f = OMEGA_H_LAMBDA(LO) {
    UnitUniformDistribution uniform_rng{seed, key, 0};
    StandardNormalDistribution normal_rng;
    for (int i = 0; i < nsamples; ++i) {
      auto uniform_x = uniform_rng();
      auto normal_x = (std::sqrt(0.2) / 2.0) * normal_rng(uniform_rng);
      auto weibull_1_x = weibull_quantile(0.5, 1.0, uniform_rng());
      auto weibull_2_x = weibull_quantile(1.0, 1.0, uniform_rng());
      auto weibull_3_x = weibull_quantile(1.5, 1.0, uniform_rng());
      auto weibull_4_x = weibull_quantile(5.0, 1.0, uniform_rng());
      contrib(d_uniform_buckets, uniform_x);
      contrib(d_normal_buckets, normal_x);
      contrib(d_weibull_1_buckets, weibull_1_x);
      contrib(d_weibull_2_buckets, weibull_2_x);
      contrib(d_weibull_3_buckets, weibull_3_x);
      contrib(d_weibull_4_buckets, weibull_4_x);
    }
  };
  parallel_for(1, f);
  auto h_uniform_buckets = HostRead<Real>(Reals(d_uniform_buckets));
  auto h_normal_buckets = HostRead<Real>(Reals(d_normal_buckets));
  auto h_weibull_1_buckets = HostRead<Real>(Reals(d_weibull_1_buckets));
  auto h_weibull_2_buckets = HostRead<Real>(Reals(d_weibull_2_buckets));
  auto h_weibull_3_buckets = HostRead<Real>(Reals(d_weibull_3_buckets));
  auto h_weibull_4_buckets = HostRead<Real>(Reals(d_weibull_4_buckets));
  // for (LO i = 0; i < nbuckets; ++i) {
  // std::cout << i << ", ";
  // std::cout << h_uniform_buckets[i] << ", ";
  // std::cout << h_normal_buckets[i] << ", ";
  // std::cout << h_weibull_1_buckets[i] << ", ";
  // std::cout << h_weibull_2_buckets[i] << ", ";
  // std::cout << h_weibull_3_buckets[i] << ", ";
  // std::cout << h_weibull_4_buckets[i];
  // std::cout << "\n";
  //}
  auto h_uniform_cdfs_w = HostWrite<Real>(nbuckets + 1);
  auto h_normal_cdfs_w = HostWrite<Real>(nbuckets + 1);
  auto h_weibull_1_cdfs_w = HostWrite<Real>(nbuckets + 1);
  auto h_weibull_2_cdfs_w = HostWrite<Real>(nbuckets + 1);
  auto h_weibull_3_cdfs_w = HostWrite<Real>(nbuckets + 1);
  auto h_weibull_4_cdfs_w = HostWrite<Real>(nbuckets + 1);
  for (LO i = 0; i < nbuckets + 1; ++i) {
    auto x = Real(i) / Real(nbuckets);
    h_uniform_cdfs_w[i] = x;
    h_normal_cdfs_w[i] =
        cumulative_general_normal_density(0.0, (std::sqrt(0.2) / 2.0), x);
    h_weibull_1_cdfs_w[i] = cumulative_weibull_density(0.5, 1.0, x);
    h_weibull_2_cdfs_w[i] = cumulative_weibull_density(1.0, 1.0, x);
    h_weibull_3_cdfs_w[i] = cumulative_weibull_density(1.5, 1.0, x);
    h_weibull_4_cdfs_w[i] = cumulative_weibull_density(5.0, 1.0, x);
  }
  // TODO: make this conversion process less terrible
  auto h_uniform_cdfs = HostRead<Real>(Reals(h_uniform_cdfs_w.write()));
  auto h_normal_cdfs = HostRead<Real>(Reals(h_normal_cdfs_w.write()));
  auto h_weibull_1_cdfs = HostRead<Real>(Reals(h_weibull_1_cdfs_w.write()));
  auto h_weibull_2_cdfs = HostRead<Real>(Reals(h_weibull_2_cdfs_w.write()));
  auto h_weibull_3_cdfs = HostRead<Real>(Reals(h_weibull_3_cdfs_w.write()));
  auto h_weibull_4_cdfs = HostRead<Real>(Reals(h_weibull_4_cdfs_w.write()));
  OMEGA_H_CHECK(std::abs(cumulative_chi_squared_density(2, 1) - 0.4) < 0.01);
  OMEGA_H_CHECK(std::abs(cumulative_chi_squared_density(3, 1) - 0.2) < 0.01);
  OMEGA_H_CHECK(std::abs(cumulative_chi_squared_density(3, 3) - 0.6) < 0.01);
  test_chi_squared(
      "uniform", h_uniform_buckets, h_uniform_cdfs, 1.0, nsamples, 0.01);
  test_chi_squared(
      "normal", h_normal_buckets, h_normal_cdfs, 2.0, nsamples, 0.01);
  test_chi_squared("Weibull(0.5, 1.0)", h_weibull_1_buckets, h_weibull_1_cdfs,
      2.0, nsamples, 0.01);
  test_chi_squared("Weibull(1.0, 1.0)", h_weibull_2_buckets, h_weibull_2_cdfs,
      2.0, nsamples, 0.01);
  test_chi_squared("Weibull(1.5, 1.0)", h_weibull_3_buckets, h_weibull_3_cdfs,
      2.0, nsamples, 0.01);
  test_chi_squared("Weibull(5.0, 1.0)", h_weibull_4_buckets, h_weibull_4_cdfs,
      2.0, nsamples, 0.01);
}
