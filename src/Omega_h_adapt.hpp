#ifndef OMEGA_H_ADAPT_HPP
#define OMEGA_H_ADAPT_HPP

#include <map>

#include <Omega_h_config.h>
#include <Omega_h_compare.hpp>
#include <Omega_h_defines.hpp>
#include <Omega_h_mark.hpp>

namespace Omega_h {

class Mesh;

struct UserTransfer {
  virtual ~UserTransfer() = default;
  virtual void out_of_line_virtual_method();
  virtual void refine(Mesh& old_mesh, Mesh& new_mesh, LOs keys2edges,
      LOs keys2midverts, Int prod_dim, LOs keys2prods, LOs prods2new_ents,
      LOs same_ents2old_ents, LOs same_ents2new_ents) = 0;
  virtual void coarsen(Mesh& old_mesh, Mesh& new_mesh, LOs keys2verts,
      Adj keys2doms, Int prod_dim, LOs prods2new_ents, LOs same_ents2old_ents,
      LOs same_ents2new_ents) = 0;
  virtual void swap(Mesh& old_mesh, Mesh& new_mesh, Int prod_dim,
      LOs keys2edges, LOs keys2prods, LOs prods2new_ents,
      LOs same_ents2old_ents, LOs same_ents2new_ents) = 0;
  virtual void swap_copy_verts(Mesh& old_mesh, Mesh& new_mesh) = 0;
};

struct TransferOpts {
  TransferOpts();
  std::map<std::string, Omega_h_Transfer> type_map;  // "density" -> CONSERVE
  std::map<std::string, std::string> integral_map;   // "density" -> "mass"
  std::map<std::string, std::string>
      velocity_density_map;  // "velocity" -> "density"
  std::map<std::string, std::string>
      velocity_momentum_map;  // "velocity" -> "momentum"
  std::map<std::string, VarCompareOpts>
      integral_diffuse_map;  // "mass" -> tolerance
  std::shared_ptr<UserTransfer> user_xfer;
  void validate(Mesh* mesh) const;
};

enum Verbosity { SILENT, EACH_ADAPT, EACH_REBUILD, EXTRA_STATS };

#ifdef OMEGA_H_USE_EGADS
struct Egads;
#endif

struct AdaptOpts {
  AdaptOpts() = default;
  AdaptOpts(Int dim);     // sets defaults
  AdaptOpts(Mesh* mesh);  // calls above
  Real min_length_desired;
  Real max_length_desired;
  Real max_length_allowed;
  Real min_quality_allowed;
  Real min_quality_desired;
  Int nsliver_layers;
  Verbosity verbosity;
  Real length_histogram_min;
  Real length_histogram_max;
  Int nlength_histogram_bins;
  Int nquality_histogram_bins;
#ifdef OMEGA_H_USE_EGADS
  Egads* egads_model;
  bool should_smooth_snap;
  Real snap_smooth_tolerance;
  bool allow_snap_failure;
#endif
  bool should_refine;
  bool should_coarsen;
  bool should_swap;
  bool should_coarsen_slivers;
  bool should_prevent_coarsen_flip;
  TransferOpts xfer_opts;
};

Real min_fixable_quality(Mesh* mesh, AdaptOpts const& opts);

/* returns false if the mesh was not modified. */
bool adapt(Mesh* mesh, AdaptOpts const& opts);

bool print_adapt_status(Mesh* mesh, AdaptOpts const& opts);
void print_adapt_histograms(Mesh* mesh, AdaptOpts const& opts);

void fix_momentum_velocity_verts(
    Mesh* mesh, std::vector<ClassPair> const& class_pairs, Int comp);

bool warp_to_limit(Mesh* mesh, AdaptOpts const& opts,
    bool exit_on_stall = false, Int max_niters = 40);
bool approach_metric(Mesh* mesh, AdaptOpts const& opts, Real min_step = 1e-4);

struct MetricSource {
  Omega_h_Source type;
  Real knob;
  std::string tag_name;
  Omega_h_Isotropy isotropy;
  Omega_h_Scales scales;
  MetricSource() = default;
  MetricSource(Omega_h_Source type_, Real knob_ = 1.0,
      std::string const& tag_name_ = "",
      Omega_h_Isotropy isotropy_ = OMEGA_H_ANISOTROPIC,
      Omega_h_Scales scales_ = OMEGA_H_SCALES);
};

struct MetricInput {
  MetricInput();
  bool verbose;
  std::vector<MetricSource> sources;
  bool should_limit_lengths;
  Real max_length;
  Real min_length;
  bool should_limit_gradation;
  Real max_gradation_rate;
  Real gradation_convergence_tolerance;
  bool should_limit_element_count;
  Real max_element_count;
  Real min_element_count;
  Real element_count_over_relaxation;
  Int nsmoothing_steps;
  void add_source(MetricSource const& src);
};

Reals generate_metrics(Mesh* mesh, MetricInput const& input);
void add_metric_tag(
    Mesh* mesh, Reals metrics, std::string const& name = "metric");
void generate_metric_tag(Mesh* mesh, MetricInput const& input);
void generate_target_metric_tag(Mesh* mesh, MetricInput const& input);
void add_implied_metric_tag(Mesh* mesh);
void add_implied_isos_tag(Mesh* mesh);
void add_implied_metric_based_on_target(Mesh* mesh);

void fix(Mesh* mesh, AdaptOpts const& adapt_opts, Omega_h_Isotropy isotropy,
    bool verbose);
void fix_for_given_metric(
    Mesh* mesh, AdaptOpts const& adapt_opts, bool verbose);

void grade_fix_adapt(
    Mesh* mesh, AdaptOpts const& opts, Reals target_metric, bool verbose);

}  // namespace Omega_h

#endif
