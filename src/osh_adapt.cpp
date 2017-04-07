#include "Omega_h_teuchos.hpp"
#include "Omega_h_cmdline.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

using namespace Omega_h;

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  auto comm = lib.world();
  CmdLine cmdline;
#ifdef OMEGA_H_USE_YAML
  auto configpath_placeholder = "input.{xml,yaml}";
#else
  auto configpath_placeholder = "input.xml";
#endif
  cmdline.add_arg<std::string>(configpath_placeholder);
  if (!cmdline.parse_all_or_help(comm, &argc, argv)) {
    return -1;
  }
  auto configpath = cmdline.get<std::string>(configpath_placeholder);
  auto pl_rcp = Teuchos::createParameterList("Omega_h");
  auto comm_teuchos = make_teuchos_comm(comm);
  update_parameters_from_file(configpath, pl_rcp.get(), comm_teuchos);
  auto inputpath = pl_rcp->get<std::string>("Input File");
  auto outputpath = pl_rcp->get<std::string>("Output File");
  auto mesh = Mesh(&lib);
  binary::read(inputpath, comm, &mesh);
  if (pl_rcp->isSublist("Target Metric")) {
    auto& target_metric_pl = pl_rcp->sublist("Target Metric");
    auto target_metric_input = MetricInput();
    update_metric_input(&target_metric_input, target_metric_pl);
    generate_target_metric_tag(&mesh, target_metric_input);
  }
  auto& metric_pl = pl_rcp->sublist("Metric");
  auto metric_input = MetricInput();
  update_metric_input(&metric_input, metric_pl);
  generate_metric_tag(&mesh, metric_input);
  auto adapt_opts = AdaptOpts(&mesh);
  if (pl_rcp->isSublist("Adapt")) {
    auto& adapt_pl = pl_rcp->sublist("Adapt");
    update_adapt_opts(&adapt_opts, adapt_pl);
  }
  do { adapt(&mesh, adapt_opts); } while (approach_metric(&mesh, adapt_opts));
  binary::write(outputpath, &mesh);
}
