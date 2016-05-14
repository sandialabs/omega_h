// poly.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
//

//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
// solve cubic equation x^3 + a * x^2 + b * x + c = 0
INLINE UInt solve_cubic(
    Real a, Real b, Real c,
    Real x[]) {
  Real a2 = square(a);
  Real q  = (a2 - 3 * b) / 9;
  Real r  = (a * (2. * a2 - 9. * b) + 27. * c) / 54;
  Real r2 = square(r);
  Real q3 = cube(q);
  Real A, B;
  if(r2 < q3) {
    Real t = r / sqrt(q3);
    t = max2(-1., t);
    t = min2( 1., t);
    t = acos(t);
    a /= 3.;
    q = -2. * sqrt(q);
    x[0] = q * cos(t / 3.) - a;
    x[1] = q * cos((t + 2. * PI) / 3.) - a;
    x[2] = q * cos((t - 2. * PI) / 3.) - a;
    return 3;
  } else {
    A = -pow(fabs(r) + sqrt(r2 - q3), 1. / 3.);
    A = fabs(A);
    B = (A == 0.) ? 0 : (q / A);
    a /= 3.;
    x[0] = (A + B) - a;
    x[1] = -0.5 * (A + B) - a;
    x[2] =  0.5 * sqrt(3.) * (A - B);
    if(fabs(x[2]) < EPSILON) {
      x[2] = x[1];
      return 2;
    }
    return 1;
  }
}
