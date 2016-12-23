#ifndef OMEGA_H_COMPARE_HPP
#define OMEGA_H_COMPARE_HPP

#include "Omega_h.hpp"
#include "Omega_h_cmdline.hpp"

namespace Omega_h {

struct VarCompareOpts {
  enum { NONE, RELATIVE, ABSOLUTE } kind;
  Real tolerance;
  Real floor;
  static VarCompareOpts zero_tolerance();
  static VarCompareOpts defaults();
  static VarCompareOpts none();
};

struct MeshCompareOpts {
  std::map<std::string, VarCompareOpts> tags2opts[4];
  VarCompareOpts time_step_opts;
  VarCompareOpts tag_opts(int dim, std::string const& name) const;
  static MeshCompareOpts init(Mesh const* mesh, VarCompareOpts var_opts);
};


Real get_real_diff(Real a, Real b, VarCompareOpts opts);
bool compare_real(Real a, Real b, VarCompareOpts opts);

Omega_h_Comparison compare_meshes(
    Mesh* a, Mesh* b, MeshCompareOpts const& opts,
    bool verbose, bool full = true);

bool check_same(Mesh* a, Mesh* b);

bool check_regression(
    std::string const& prefix, Mesh* mesh);

void get_diff_program_cmdline(
    std::string const& a_name,
    std::string const& b_name,
    CmdLine* p_cmdline);

void accept_diff_program_cmdline(CmdLine const& cmdline,
    Mesh const* mesh,
    MeshCompareOpts* p_opts,
    Omega_h_Comparison* p_max_result);

}

#endif
