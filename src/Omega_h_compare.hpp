#ifndef OMEGA_H_COMPARE_HPP
#define OMEGA_H_COMPARE_HPP

#include "Omega_h.hpp"
#include "Omega_h_cmdline.hpp"

namespace Omega_h {

struct ArrayCompareOpts {
  enum { RELATIVE, ABSOLUTE } kind;
  Real tolerance;
  Real floor;
};

struct MeshCompareOpts {
  std::map<std::string, ArrayCompareOpts> tags2opts[4];
  ArrayCompareOpts default_tag_opts;
};

ArrayCompareOpts get_tag_opts(MeshCompareOpts const& opts,
    int dim, std::string const& name);

MeshCompareOpts get_exodiff_defaults();
MeshCompareOpts get_zero_tolerance();


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
    MeshCompareOpts* p_opts,
    Omega_h_Comparison* p_max_result);

}

#endif
