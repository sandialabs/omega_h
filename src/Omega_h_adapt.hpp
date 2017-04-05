#ifndef OMEGA_H_ADAPT_HPP
#define OMEGA_H_ADAPT_HPP

#include <map>

#include <Omega_h_config.h>
#include <Omega_h_compare.hpp>
#include <Omega_h_defines.hpp>

namespace Omega_h {

class Mesh;

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
  bool should_conserve_size;
  void validate(Mesh* mesh) const;
};

enum Verbosity { SILENT, EACH_ADAPT, EACH_REBUILD, EXTRA_STATS };

#ifdef OMEGA_H_USE_EGADS
struct Egads;
#endif

struct AdaptOpts {
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
#endif
  Int max_motion_steps;
  Real motion_step_size;
  bool should_refine;
  bool should_coarsen;
  bool should_swap;
  bool should_coarsen_slivers;
  bool should_move_for_quality;
  bool should_allow_pinching;
  TransferOpts xfer_opts;
};

Real min_fixable_quality(Mesh* mesh, AdaptOpts const& opts);

/* returns false if the mesh was not modified. */
bool adapt(Mesh* mesh, AdaptOpts const& opts);

bool print_adapt_status(Mesh* mesh, AdaptOpts const& opts);
void print_adapt_histograms(Mesh* mesh, AdaptOpts const& opts);

void fix_momentum_velocity_verts(
    Mesh* mesh, Int class_dim, I32 class_id, Int comp);

bool warp_to_limit(Mesh* mesh, AdaptOpts const& opts);
bool approach_metric(Mesh* mesh, AdaptOpts const& opts);

struct MetricSource {
  Omega_h_Source kind;
  bool should_scale;
  std::string tag_name;
  Real knob;
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
};

Reals generate_metrics(Mesh* mesh, MetricInput const& input);
void add_metric_tag(
    Mesh* mesh, Reals metrics, std::string const& name = "metric");
void generate_metric_tag(Mesh* mesh, MetricInput const& input);
void generate_target_metric_tag(Mesh* mesh, MetricInput const& input);
void add_implied_metric_tag(Mesh* mesh);

}  // namespace Omega_h

#endif
