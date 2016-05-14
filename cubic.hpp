// poly.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
//

// solve cubic equation roots^3 + a * roots^2 + b * roots + c = 0
INLINE UInt solve_cubic(
    Real a, Real b, Real c,
    Real roots[], UInt mults[]) {
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
    roots[0] = q * cos(t / 3.) - a;
    roots[1] = q * cos((t + 2. * PI) / 3.) - a;
    roots[2] = q * cos((t - 2. * PI) / 3.) - a;
    mults[0] = 1;
    mults[1] = 1;
    mults[2] = 1;
    return 3;
  } else {
    A = -pow(fabs(r) + sqrt(r2 - q3), 1. / 3.);
    A = fabs(A);
    B = (A == 0.) ? 0. : (q / A);
    a /= 3.;
    roots[0] = (A + B) - a;
    Real im =  0.5 * sqrt(3.) * (A - B);
    if (fabs(im) < EPSILON) {
      if (are_close(roots[0], roots[1])) {
        mults[0] = 3;
        return 1;
      }
      mults[0] = 1;
      mults[1] = 2;
      roots[1] = -0.5 * (A + B) - a;
      return 2;
    }
    return 1;
  }
}
