#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_vtk.hpp"
#include "Omega_h_xml.hpp"

#include <iostream>
#include <sstream>

using namespace Omega_h;

static void test_file_components() {
  using namespace binary;
  std::stringstream stream;
  std::string s = "foo";
  LO n = 10;
#ifdef OMEGA_H_USE_ZLIB
  bool is_compressed = true;
#else
  bool is_compressed = false;
#endif
  I8 a = 2;
  write_value(stream, a);
  I32 b = 42 * 1000 * 1000;
  write_value(stream, b);
  I64 c = I64(42) * 1000 * 1000 * 1000;
  write_value(stream, c);
  Real d = 4.2;
  write_value(stream, d);
  Read<I8> aa(n, 0, a);
  write_array(stream, aa);
  Read<I32> ab(n, 0, b);
  write_array(stream, ab);
  Read<I64> ac(n, 0, c);
  write_array(stream, ac);
  Read<Real> ad(n, 0, d);
  write_array(stream, ad);
  write(stream, s);
  I8 a2;
  read_value(stream, a2);
  OMEGA_H_CHECK(a == a2);
  I32 b2;
  read_value(stream, b2);
  OMEGA_H_CHECK(b == b2);
  I64 c2;
  read_value(stream, c2);
  OMEGA_H_CHECK(c == c2);
  Real d2;
  read_value(stream, d2);
  OMEGA_H_CHECK(d == d2);
  Read<I8> aa2;
  read_array(stream, aa2, is_compressed);
  OMEGA_H_CHECK(aa2 == aa);
  Read<I32> ab2;
  read_array(stream, ab2, is_compressed);
  OMEGA_H_CHECK(ab2 == ab);
  Read<I64> ac2;
  read_array(stream, ac2, is_compressed);
  OMEGA_H_CHECK(ac2 == ac);
  Read<Real> ad2;
  read_array(stream, ad2, is_compressed);
  OMEGA_H_CHECK(ad2 == ad);
  std::string s2;
  read(stream, s2);
  OMEGA_H_CHECK(s == s2);
}

static void build_empty_mesh(Mesh* mesh, Int dim) {
  build_from_elems_and_coords(mesh, OMEGA_H_SIMPLEX, dim, LOs({}), Reals({}));
}

static void test_file(Library* lib, Mesh* mesh0) {
  std::stringstream stream;
  binary::write(stream, mesh0);
  Mesh mesh1(lib);
  mesh1.set_comm(lib->self());
  binary::read(stream, &mesh1, binary::latest_version);
  mesh1.set_comm(lib->world());
  auto opts = MeshCompareOpts::init(mesh0, VarCompareOpts::zero_tolerance());
  compare_meshes(mesh0, &mesh1, opts, true, true);
  OMEGA_H_CHECK(*mesh0 == mesh1);
}

static void test_file(Library* lib) {
  {
    auto mesh0 = build_box(lib->world(), OMEGA_H_SIMPLEX, 1., 1., 1., 1, 1, 1);
    test_file(lib, &mesh0);
  }
  {
    Mesh mesh0(lib);
    build_empty_mesh(&mesh0, 3);
    test_file(lib, &mesh0);
  }
}

static void test_xml() {
  xml::Tag tag;
  OMEGA_H_CHECK(!xml::parse_tag("AQAAAAAAAADABg", &tag));
  OMEGA_H_CHECK(!xml::parse_tag("   <Foo bar=\"qu", &tag));
  OMEGA_H_CHECK(!xml::parse_tag("   <Foo bar=", &tag));
  OMEGA_H_CHECK(xml::parse_tag("   <Foo bar=\"quux\"   >", &tag));
  OMEGA_H_CHECK(tag.elem_name == "Foo");
  OMEGA_H_CHECK(tag.attribs["bar"] == "quux");
  OMEGA_H_CHECK(tag.type == xml::Tag::START);
  OMEGA_H_CHECK(xml::parse_tag("   <Elem att=\"val\"  answer=\"42\" />", &tag));
  OMEGA_H_CHECK(tag.elem_name == "Elem");
  OMEGA_H_CHECK(tag.attribs["att"] == "val");
  OMEGA_H_CHECK(tag.attribs["answer"] == "42");
  OMEGA_H_CHECK(tag.type == xml::Tag::SELF_CLOSING);
  OMEGA_H_CHECK(xml::parse_tag("</Foo>", &tag));
  OMEGA_H_CHECK(tag.elem_name == "Foo");
  OMEGA_H_CHECK(tag.type == xml::Tag::END);
}

static void test_read_vtu(Mesh* mesh0) {
  std::stringstream stream;
  vtk::write_vtu(stream, mesh0, mesh0->dim(), vtk::get_all_vtk_tags(mesh0));
  Mesh mesh1(mesh0->library());
  vtk::read_vtu(stream, mesh0->comm(), &mesh1);
  auto opts = MeshCompareOpts::init(mesh0, VarCompareOpts::zero_tolerance());
  OMEGA_H_CHECK(
      OMEGA_H_SAME == compare_meshes(mesh0, &mesh1, opts, true, false));
}

static void test_read_vtu(Library* lib) {
  auto mesh0 = build_box(lib->world(), OMEGA_H_SIMPLEX, 1., 1., 1., 1, 1, 1);
  test_read_vtu(&mesh0);
}

int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv);
  OMEGA_H_CHECK(std::string(lib.version()) == OMEGA_H_SEMVER);
  test_file_components();
  test_file(&lib);
  test_xml();
  test_read_vtu(&lib);
}
