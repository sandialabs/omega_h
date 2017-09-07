#include "Omega_h_quality.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>

static void tet_run(Omega_h::Int dihedral_angle_in_degrees, std::ostream& os) {
  auto dihedral_angle = Omega_h::Real(dihedral_angle_in_degrees) / 180. * Omega_h::PI;
  Omega_h::Few<Omega_h::Vector<3>, 4> x;
  x[0] = Omega_h::vector_3( 0.5, 0, 0);
  x[2] = Omega_h::vector_3(-0.5, 0, 0);
  auto h = std::sqrt(3.0) / 2.0;
  x[1] = Omega_h::vector_3(0, h, 0);
  x[3] = Omega_h::vector_3(0, std::cos(dihedral_angle) * h, std::sin(dihedral_angle) * h);
  auto m = Omega_h::matrix_1x1(1.0);
  auto quality = Omega_h::metric_element_quality<3>(x, m);
  os << quality;
}

static void tri_run(Omega_h::Int dihedral_angle_in_degrees, std::ostream& os) {
  auto dihedral_angle = Omega_h::Real(dihedral_angle_in_degrees) / 180. * Omega_h::PI;
  auto half_dihedral_angle = dihedral_angle / 2.;
  auto half_cord_length = sin(half_dihedral_angle);
  auto cord_length = 2. * half_cord_length;
  auto height = sqrt(1. - Omega_h::square(half_cord_length));
  auto area = cord_length * height / 2.;
  auto msl = (Omega_h::square(cord_length) + 2.) / 3.;
  auto quality = Omega_h::mean_ratio<2>(area, msl);
  os << quality;
}

int main() {
  {
  std::ofstream plotfile("dihedral_quality.csv");
  plotfile << std::scientific << std::setprecision(16);
  for (Omega_h::Int angle = 80; angle >= 0; --angle) {
    plotfile << angle << ", ";
    tri_run(angle, plotfile);
    plotfile << ", ";
    tet_run(angle, plotfile);
    plotfile << '\n';
  }
  }
//constexpr Int ndihedral_samples = 72;
//{
//std::ofstream hppfile("Omega_h_dihedral_quality.cpp");
//hppfile << "#ifndef OMEGA_H_DIHEDRAL_QUALITY_HPP\n";
//hppfile << "#define OMEGA_H_DIHEDRAL_QUALITY_HPP\n\n";
//hppfile << "#include <Omega_h_defines.hpp>\n\n";
//hppfile << "namespace Omega_h {\n\n";
//hppfile << "static constexpr Int ndihedral_samples = " << ndihedral_samples << ";\n";
//hppfile << "extern const Real best_tri_qual_for_dihedral[ndihedral_samples];\n";
//hppfile << "extern const Real best_tet_qual_for_dihedral[ndihedral_samples];\n\n";
//hppfile << "}  // end namespace Omega_h\n\n";
//hppfile << "#endif\n";
//}
//{
//std::ofstream cppfile("Omega_h_dihedral_quality.cpp");
//cppfile << std::scientific << std::setprecision(16);
//cppfile << "#include \"Omega_h_dihedral_quality.hpp\"\n\n";
//cppfile << "namespace Omega_h {\n\n";
//cppfile << "const Real best_tri_qual_for_dihedral[ndihedral_samples] = {\n";
//for (Omega_h::Int angle = 0; angle < ndihedral_samples; ++angle) {
//  tri_run(angle, cppfile);
//  if (angle + 1 < ndihedral_samples) cppfile << ',';
//  cppfile << '\n';
//}
//cppfile << "};\n\n";
//cppfile << "const Real best_tri_qual_for_dihedral[ndihedral_samples] = {\n";
//}
}
