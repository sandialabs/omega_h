// poly.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
//

// solve cubic equation roots^3 + a * roots^2 + b * roots + c = 0
INLINE UInt solve_cubic(
    Real a, Real b, Real c,
    Real roots[], UInt mults[]) {
// https://en.wikipedia.org/wiki/Cubic_function#Reduction_to_a_depressed_cubic
  // this is b^2 from wikipedia
  Real a2 = square(a);
  // this is (-p / 3) from wikipedia
  Real q  = (a2 - 3 * b) / 9;
  // this is (q/2) from wikipedia
  Real r  = (a * (2. * a2 - 9. * b) + 27. * c) / 54;
  // this is (q^2/4) from wikipedia
  Real r2 = square(r);
  // this is (-p^3/27) from wikipedia
  Real q3 = cube(q);
  //this is the (opposite of the) Cardano assumption ((q^2/4) + (p^3/27)) > 0
  if(r2 < q3) {
// https://en.wikipedia.org/wiki/Cubic_function#Three_real_roots
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
// https://en.wikipedia.org/wiki/Cubic_function#Cardano.27s_method
    Real u3 = -r - sqrt(r2 - q3);
    std::cerr << "u3 " << u3 << '\n';
    Real u = pow(u3, 1./3.);
    std::cerr << "u " << u << '\n';
    Real t1 = u + v;
    // recall x = t - (b/(3a)), in our case x = t - (a/3)
    roots[0] = t1 - (a / 3.);
    std::cerr << "roots[0] " << roots[0] << '\n';
    Real t_real = -0.5 * (u + v);
    Real t_imag = 0.5 * sqrt(3) * (u - v);
    std::cerr << "t_imag " << t_imag << '\n';
    roots[1] = (t_real) - (a / 3.);
    if (fabs(t_imag) < EPSILON) {
      std::cerr << "roots[1] " << roots[1] << '\n';
      if (are_close(roots[0], roots[1])) {
        mults[0] = 3;
        return 1;
      }
      mults[0] = 1;
      mults[1] = 2;
      return 2;
    }
    mults[0] = 1;
    return 1;
  }
}
