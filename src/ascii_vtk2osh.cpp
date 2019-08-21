#include <fstream>
#include "Omega_h_build.hpp"
#include "Omega_h_cmdline.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_element.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_mesh.hpp"

// this is super-specific to a single kind of output generated
// by tetgen from CT data with no attempt for generalization.

static std::string get_line(std::istream& stream) {
  std::string line;
  std::getline(stream, line);
  return line;
}

static void check_line(std::istream& stream, std::string want) {
  auto line = get_line(stream);
  if (line != want) {
    Omega_h_fail("wanted line: \"%s\" but got: \"%s\"\n",
        line.c_str(), want.c_str());
  }
  OMEGA_H_CHECK(line == want);
}

static void check_header(std::istream& stream) {
  check_line(stream, "# vtk DataFile Version 2.0");
  check_line(stream, "Unstructured Grid");
  check_line(stream, "ASCII");
  check_line(stream, "DATASET UNSTRUCTURED_GRID");
}

template <class T>
static void read_val(std::istream& stream, T& value) {
  stream >> value;
}

static void eat_newlines(std::istream& stream) {
  while (stream.peek() == int('\n')) stream.get();
}

static int get_npoints(std::istream& stream) {
  int npoints;
  std::string points_marker;
  std::string double_marker;
  read_val(stream, points_marker);
  read_val(stream, npoints);
  read_val(stream, double_marker);
  OMEGA_H_CHECK(points_marker == "POINTS");
  OMEGA_H_CHECK(double_marker == "double");
  return npoints;
}

static Omega_h::HostWrite<Omega_h::Real> get_coords(std::istream& stream) {
  int npoints = get_npoints(stream);
  Omega_h::HostWrite<Omega_h::Real> coords(npoints * 3);
  for (int pt = 0; pt < npoints; ++pt) {
    Omega_h::Real val;
    for (int d = 0; d < 3; ++d) {
      read_val(stream, val);
      coords[pt * 3 + d] = val;
    }
  }
  eat_newlines(stream);
  return coords;
}

static int get_cell_nelems(std::istream& stream) {
  int nelems;
  int garbage;
  std::string cells_marker;
  read_val(stream, cells_marker);
  read_val(stream, nelems);
  read_val(stream, garbage);
  OMEGA_H_CHECK(cells_marker == "CELLS");
  OMEGA_H_CHECK(nelems > 0);
  return nelems;
}

static Omega_h::HostWrite<Omega_h::LO> get_ev2v(std::istream& stream) {
  int nelems = get_cell_nelems(stream);
  int neev = Omega_h::element_degree(OMEGA_H_SIMPLEX,  3, 0);
  Omega_h::HostWrite<Omega_h::LO> ev2v(nelems * neev);
  for (int elem = 0; elem < nelems; ++elem) {
    Omega_h::LO val;
    read_val(stream, val);
    OMEGA_H_CHECK(val == neev);
    for (int v = 0; v < neev; ++v) {
      read_val(stream, val);
      ev2v[elem * neev + v] = val;
    }
  }
  eat_newlines(stream);
  return ev2v;
}

static void eat_elem_types(std::istream& stream) {
  int nlines;
  int garbage;
  std::string cell_type_marker;
  read_val(stream, cell_type_marker);
  read_val(stream, nlines);
  OMEGA_H_CHECK(cell_type_marker == "CELL_TYPES");
  for (int i = 0; i < nlines; ++i) read_val(stream, garbage);
  eat_newlines(stream);
}

static int get_mat_id_nelems(std::istream& stream) {
  int nelems;
  std::string cell_data_marker;
  std::string scalars_marker;
  read_val(stream, cell_data_marker);
  read_val(stream, nelems);
  read_val(stream, scalars_marker);
  OMEGA_H_CHECK(cell_data_marker == "CELL_DATA");
  get_line(stream);
  get_line(stream);
  return nelems;
}

static Omega_h::HostWrite<Omega_h::LO> get_elem_mat_ids(std::istream& stream) {
  int nelems = get_mat_id_nelems(stream);
  Omega_h::HostWrite<Omega_h::LO> mat(nelems);
  for (int elem = 0; elem < nelems; ++elem) {
    Omega_h::LO mat_id;
    read_val(stream, mat_id);
    mat[elem] = mat_id;
  }
  eat_newlines(stream);
  return mat;
}

static void build_mesh(
    Omega_h::Mesh* mesh,
    Omega_h::HostWrite<Omega_h::Real> h_coords,
    Omega_h::HostWrite<Omega_h::LO> h_ev2v) {
  Omega_h::Read<Omega_h::Real> coords(h_coords.write());
  Omega_h::Read<Omega_h::LO> ev2v(h_ev2v.write());
  Omega_h::build_from_elems_and_coords(
      mesh, OMEGA_H_SIMPLEX, 3, ev2v, coords);
}

static void classify(
    Omega_h::Mesh* mesh,
    Omega_h::HostWrite<Omega_h::LO> h_mat) {
  Omega_h::Read<Omega_h::LO> mat(h_mat.write());
  mesh->add_tag(3, "class_id", 1, mat);
  Omega_h::finalize_classification(mesh);
}

static void build(Omega_h::Mesh* mesh, std::string vtk_path) {
  std::ifstream file(vtk_path.c_str());
  OMEGA_H_CHECK(file.is_open());
  check_header(file);
  auto coords = get_coords(file);
  auto ev2v = get_ev2v(file);
  eat_elem_types(file);
  auto mat = get_elem_mat_ids(file);
  OMEGA_H_CHECK(file.eof());
  build_mesh(mesh, coords, ev2v);
  classify(mesh, mat);
  Omega_h::reorder_by_hilbert(mesh);
}

int main(int argc, char** argv) {
  OMEGA_H_CHECK(argc == 3);
  auto lib = Omega_h::Library(&argc, &argv);
  auto comm = lib.world();
  auto cmdline = Omega_h::CmdLine();
  cmdline.add_arg<std::string>("mesh.vtk");
  cmdline.add_arg<std::string>("out.osh");
  if (!cmdline.parse_final(comm, &argc, argv)) return -1;
  auto vtk_path = cmdline.get<std::string>("mesh.vtk");
  auto osh_path = cmdline.get<std::string>("out.osh");
  Omega_h::Mesh mesh(&lib);
  build(&mesh, vtk_path);
  Omega_h::binary::write(osh_path, &mesh);
}
