#include <Omega_h_adapt.hpp>
#include <Omega_h_cmdline.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>

#ifdef OMEGA_H_USE_EGADS
#include <Omega_h_egads.hpp>
#endif

#include <iostream>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::CmdLine cmdline;
  cmdline.add_arg<std::string>("mesh_in.osh");
  cmdline.add_arg<std::string>("mesh_out.osh");
#ifdef OMEGA_H_USE_EGADS
  auto& model_flag = cmdline.add_flag("--model", "optional EGADS model");
  model_flag.add_arg<std::string>("model.step");
#endif
  cmdline.add_flag("--isotropic", "aim for an isotropic mesh");
  if (!cmdline.parse_final(lib.world(), &argc, argv)) return -1;
  auto path_in = cmdline.get<std::string>("mesh_in.osh");
  auto path_out = cmdline.get<std::string>("mesh_out.osh");
  Omega_h::Mesh mesh(&lib);
  std::cout << "reading in " << path_in << '\n';
  Omega_h::binary::read(path_in, lib.world(), &mesh);
#ifdef OMEGA_H_USE_EGADS
  Omega_h::Egads* eg = nullptr;
  if (cmdline.parsed("--model")) {
    auto model_path = cmdline.get<std::string>("--model", "model.step");
    std::cout << "reading in " << model_path << '\n';
    eg = Omega_h::egads_load(model_path);
  }
#endif
  Omega_h::AdaptOpts opts(&mesh);
#ifdef OMEGA_H_USE_EGADS
  opts.egads_model = eg;
#endif
  Omega_h_Isotropy isotropy = cmdline.parsed("--isotropic") ? OMEGA_H_ISO_LENGTH : OMEGA_H_ANISOTROPIC;
  Omega_h::fix(&mesh, opts, isotropy, true);
  mesh.remove_tag(Omega_h::VERT, "metric");
  std::cout << "writing out " << path_out << '\n';
  Omega_h::binary::write(path_out, &mesh);
#ifdef OMEGA_H_USE_EGADS
  if (eg != nullptr) {
    Omega_h::egads_free(eg);
  }
#endif
}
