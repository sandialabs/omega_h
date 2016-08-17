#include "quality.hpp"

#include <iostream>

/* This program exists to determine the worst-case
 * number of tetrahedra around a vertex and its relationship
 * to the mean ratio quality measure.
 * The worst-case solid angle determines the worst-case
 * number of tetrahedra around a single vertex.
 *
 * A needle tetrahedron is formed using the center of the unit
 * sphere and an equilateral triangle on the surface of the unit sphere.
 * This should maximize quality over solid angle.
 *
 * The program takes as argument the small angle (in degrees) of one of
 * the isosceles triangles from the center to the surface.
 */

/* http://www.had2know.com/academics/dihedral-angle-calculator-polyhedron.html
 * https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
 * https://en.wikipedia.org/wiki/Tetrahedron#Volume
 */

static void tet_run(osh::Real side_angle_in_degrees) {
  auto side_angle = side_angle_in_degrees / 180. * osh::PI;
  std::cout << "side_angle " << side_angle << '\n';
  auto dihedral_angle = acos((cos(side_angle) - osh::square(cos(side_angle))) /
                             (osh::square(sin(side_angle))));
  std::cout << "dihedral_angle " << dihedral_angle << '\n';
  std::cout << "dihedral_angle in degrees " << (dihedral_angle * 180. / osh::PI)
            << '\n';
  auto solid_angle = 3. * dihedral_angle - osh::PI;
  std::cout << "solid_angle " << solid_angle << '\n';
  auto degree = 4. * osh::PI / solid_angle;
  std::cout << "degree " << degree << '\n';
  auto half_side_angle = side_angle / 2.;
  auto half_cord_length = sin(half_side_angle);
  auto cord_length = 2. * half_cord_length;
  std::cout << "cord_length " << cord_length << '\n';
  auto surf_tri_height = sqrt(3.) / 2. * cord_length;
  auto inradius = surf_tri_height / 3.;
  auto circumradius = inradius * 2.;
  std::cout << "circumradius " << circumradius << '\n';
  auto surf_tri_area = sqrt(3.) / 4. * osh::square(cord_length);
  auto tet_height = sqrt(1. - osh::square(circumradius));
  auto tet_volume = 1. / 3. * surf_tri_area * tet_height;
  std::cout << "tet_volume " << tet_volume << '\n';
  auto msl = (3 * osh::square(cord_length) + 3.) / 6.;
  auto quality = osh::tet_mean_ratio(tet_volume, msl);
  std::cout << "quality " << quality << '\n';
}

/* Here is the 2D triangle analogue.
 */

static void tri_run(osh::Real side_angle_in_degrees) {
  auto side_angle = side_angle_in_degrees / 180. * osh::PI;
  std::cout << "side_angle " << side_angle << '\n';
  auto degree = 2. * osh::PI / side_angle;
  std::cout << "degree " << degree << '\n';
  auto half_side_angle = side_angle / 2.;
  auto half_cord_length = sin(half_side_angle);
  auto cord_length = 2. * half_cord_length;
  std::cout << "cord_length " << cord_length << '\n';
  auto height = sqrt(1. - osh::square(half_cord_length));
  auto area = cord_length * height / 2.;
  std::cout << "area " << area << '\n';
  auto msl = (osh::square(cord_length) + 2.) / 3.;
  auto quality = osh::tri_mean_ratio(area, msl);
  std::cout << "quality " << quality << '\n';
}

int main(int argc, char** argv) {
  OMEGA_H_CHECK(argc == 2);
  std::cout << "\nTETRAHEDRON:\n";
  tet_run(atof(argv[1]));
  std::cout << "\nTRIANGLE:\n";
  tri_run(atof(argv[1]));
}
