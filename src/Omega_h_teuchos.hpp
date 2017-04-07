#ifndef OMEGA_H_TEUCHOS_HPP
#define OMEGA_H_TEUCHOS_HPP

#include <Omega_h_adapt.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp>
#if defined(HAVE_TEUCHOSCORE_YAML_CPP) && HAVE_TEUCHOSCORE_YAML_CPP
#define OMEGA_H_USE_YAML
#endif

namespace Omega_h {

template <typename T>
void set_if_given(T* var, Teuchos::ParameterList const& pl, std::string const& name) {
  if (pl.isType<T>(name)) *var = pl.get<T>(name);
}

void update_var_compare_opts(VarCompareOpts* opts, Teuchos::ParameterList const& pl);
void update_transfer_opts(TransferOpts* opts, Teuchos::ParameterList const& pl);
void update_adapt_opts(AdaptOpts* opts, Teuchos::ParameterList const& pl);
MetricSource get_metric_source(Teuchos::ParameterList const& pl);
void update_metric_input(MetricInput* input, Teuchos::ParameterList const& pl);

Teuchos::RCP<Teuchos::Comm<int>> make_teuchos_comm(CommPtr comm_osh);
void update_parameters_from_file(
    std::string const& filepath, Teuchos::RCP<Teuchos::ParameterList> pl,
    Teuchos::RCP<Teuchos::Comm<int>> comm);

void write_parameters(
    std::string const& filepath, Teuchos::ParameterList const& pl);

}

#endif
