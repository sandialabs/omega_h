#include <Omega_h_amr.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_map.hpp>

template <int dim>
OMEGA_H_INLINE double eval_rc(Omega_h::Vector<dim> c);

template <>
OMEGA_H_INLINE double eval_rc<2>(Omega_h::Vector<2> c) {
  auto rc2 = (c[0] - 0.5) * (c[0] - 0.5) + (c[1] - 0.5) * (c[1] - 0.5);
  return std::sqrt(rc2);
}

template <>
OMEGA_H_INLINE double eval_rc<3>(Omega_h::Vector<3> c) {
  auto rc2 = (c[0] - 0.5) * (c[0] - 0.5) + (c[1] - 0.5) * (c[1] - 0.5) +
             (c[2] - 0.5) * (c[2] - 0.5);
  return std::sqrt(rc2);
}

template <int dim>
static Omega_h::Bytes mark(Omega_h::Mesh* m, int level) {
  auto coords = m->coords();
  auto mids = Omega_h::average_field(m, dim, dim, coords);
  auto is_leaf = m->ask_leaves(dim);
  auto leaf_elems = Omega_h::collect_marked(is_leaf);
  Omega_h::Write<Omega_h::Byte> marks(m->nelems(), 0);
  auto f = OMEGA_H_LAMBDA(Omega_h::LO e) {
    auto elem = leaf_elems[e];
    auto c = Omega_h::get_vector<dim, Omega_h::Reals>(mids, elem);
    auto rc = eval_rc<dim>(c);
    auto r = 0.25;
    auto tol = 0.314 / static_cast<double>(level);
    if (std::abs(rc - r) < tol) marks[elem] = 1;
  };
  Omega_h::parallel_for(leaf_elems.size(), f);
  return Omega_h::enforce_one_level(m, dim - 1, marks);
}

static void run_2D_adapt(Omega_h::Library* lib) {
  auto w = lib->world();
  auto f = OMEGA_H_HYPERCUBE;
  auto m = Omega_h::build_box(w, f, 1.0, 1.0, 0.0, 2, 2, 0);
  Omega_h::vtk::Writer writer("out_amr_2D", &m);
  writer.write();
  for (int i = 1; i < 5; ++i) {
    auto xfer_opts = Omega_h::TransferOpts();
    auto marks = mark<2>(&m, i);
    Omega_h::amr_refine(&m, marks, xfer_opts);
    writer.write();
  }
}

static void run_3D_adapt(Omega_h::Library* lib) {
  auto w = lib->world();
  auto f = OMEGA_H_HYPERCUBE;
  auto m = Omega_h::build_box(w, f, 1.0, 1.0, 1.0, 2, 2, 2);
  Omega_h::vtk::Writer writer("out_amr_3D", &m);
  writer.write();
  for (int i = 1; i < 5; ++i) {
    auto xfer_opts = Omega_h::TransferOpts();
    auto marks = mark<3>(&m, i);
    Omega_h::amr_refine(&m, marks, xfer_opts);
    writer.write();
  }
}

int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  run_2D_adapt(&lib);
  run_3D_adapt(&lib);
}
