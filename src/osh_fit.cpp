#include <Omega_h_adapt.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_recover.hpp>

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  // build mesh
  auto mesh = Omega_h::build_box(
      lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 16, 16, 0);
  // attach density field
  auto coords = mesh.coords();
  auto u_w = Omega_h::Write<Omega_h::Real>(mesh.nverts());
  constexpr auto dim = 2;
  auto f = OMEGA_H_LAMBDA(Omega_h::LO v) {
    auto x = Omega_h::get_vector<dim>(coords, v);
    auto r = Omega_h::norm(x - Omega_h::fill_vector<dim>(0.5));
    u_w[v] = (r < 0.4 ? 1.0 : 0.0);
  };
  Omega_h::parallel_for(mesh.nverts(), f);
  auto u = Omega_h::Reals(u_w);
  mesh.add_tag(OMEGA_H_VERT, "u", 1, u);
  Omega_h::vtk::write_vtu("u.vtu", &mesh);
  // adapt the mesh
  auto grad_u = Omega_h::recover_gradients(&mesh, u);
  auto target_metric = Omega_h::get_gradient_metrics(dim, grad_u, 5.0e-2);
  target_metric =
      Omega_h::clamp_metrics(mesh.nverts(), target_metric, 0.0, 1.0);
  mesh.add_tag(
      OMEGA_H_VERT, "target_metric", Omega_h::symm_ncomps(dim), target_metric);
  Omega_h::vtk::write_vtu("target.vtu", &mesh);
  auto adapt_opts = Omega_h::AdaptOpts(&mesh);
  adapt_opts.xfer_opts.type_map["u"] = OMEGA_H_LINEAR_INTERP;
  Omega_h::grade_fix_adapt(&mesh, adapt_opts, target_metric, true);
  Omega_h::vtk::write_vtu("adapted.vtu", &mesh);
  // start the classification process
  u = mesh.get_array<Omega_h::Real>(OMEGA_H_VERT, "u");
  auto elem_u = Omega_h::average_field(&mesh, mesh.dim(), 1, u);
  Omega_h::ClassId inside_class_id = 42;
  Omega_h::ClassId outside_class_id = 31;
  Omega_h::ClassId interface_class_id = 8;
  auto elem_class_ids_w = Omega_h::deep_copy(
      mesh.get_array<Omega_h::ClassId>(mesh.dim(), "class_id"));
  auto f2 = OMEGA_H_LAMBDA(Omega_h::LO e) {
    elem_class_ids_w[e] =
        (elem_u[e] >= 0.5 ? inside_class_id : outside_class_id);
  };
  Omega_h::parallel_for(mesh.nelems(), f2);
  auto elem_class_ids = Omega_h::Read<Omega_h::ClassId>(elem_class_ids_w);
  mesh.set_tag(mesh.dim(), "class_id", elem_class_ids);
  auto side_class_ids_w = Omega_h::deep_copy(
      mesh.get_array<Omega_h::ClassId>(mesh.dim() - 1, "class_id"));
  auto side_class_dims =
      mesh.get_array<Omega_h::Byte>(mesh.dim() - 1, "class_dim");
  auto f3 = OMEGA_H_LAMBDA(Omega_h::LO s) {
    side_class_ids_w[s] =
        side_class_dims[s] > (dim - 1) ? -1 : side_class_ids_w[s];
  };
  Omega_h::parallel_for(mesh.nents(mesh.dim() - 1), f3);
  auto side_class_ids = Omega_h::Read<Omega_h::ClassId>(side_class_ids_w);
  mesh.set_tag(mesh.dim() - 1, "class_id", side_class_ids);
  auto relaxed_projection = false;
  Omega_h::project_classification(&mesh, mesh.dim() - 1, relaxed_projection);
  side_class_ids_w = Omega_h::deep_copy(
      mesh.get_array<Omega_h::ClassId>(mesh.dim() - 1, "class_id"));
  auto f4 = OMEGA_H_LAMBDA(Omega_h::LO s) {
    side_class_ids_w[s] =
        side_class_ids_w[s] == -1 ? interface_class_id : side_class_ids_w[s];
  };
  Omega_h::parallel_for(mesh.nents(mesh.dim() - 1), f4);
  side_class_ids = Omega_h::Read<Omega_h::ClassId>(side_class_ids_w);
  mesh.set_tag(mesh.dim() - 1, "class_id", side_class_ids);
  auto vert_class_ids_w = Omega_h::deep_copy(
      mesh.get_array<Omega_h::ClassId>(OMEGA_H_VERT, "class_id"));
  auto vert_class_dims =
      mesh.get_array<Omega_h::Byte>(OMEGA_H_VERT, "class_dim");
  auto f5 = OMEGA_H_LAMBDA(Omega_h::LO v) {
    vert_class_ids_w[v] = vert_class_dims[v] > 0 ? -1 : vert_class_ids_w[v];
  };
  Omega_h::parallel_for(mesh.nverts(), f5);
  auto vert_class_ids = Omega_h::Read<Omega_h::ClassId>(vert_class_ids_w);
  mesh.set_tag(OMEGA_H_VERT, "class_id", vert_class_ids);
  Omega_h::project_classification(&mesh, OMEGA_H_VERT, relaxed_projection);
  Omega_h::vtk::write_vtu("reclassified.vtu", &mesh);
  Omega_h::vtk::write_vtu("reclassified_edges.vtu", &mesh, 1);
}
