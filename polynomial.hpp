// poly.cpp : solution of cubic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
//
// Dan Ibanez: took initial code, fixed lots of corner cases,
// added geometric multiplicities output, focus only on real roots

// solve cubic equation x^3 + a * x^2 + b * x + c = 0
INLINE UInt solve_cubic(
    Real a, Real b, Real c,
    Few<Real, 3>& roots, Few<UInt, 3>& mults) {
// http://mathworld.wolfram.com/CubicFormula.html
// https://en.wikipedia.org/wiki/Cubic_function#Reduction_to_a_depressed_cubic
  // this is b^2 from wikipedia
  Real a2 = square(a);
  // this is (-p/3) from wikipedia, Q from wolfram mathworld
  Real q  = (a2 - 3. * b) / 9.;
  // this is (q/2) from wikipedia, -R from wolfram mathworld
  Real r  = (a * (2. * a2 - 9. * b) + 27. * c) / 54.;
  // this is (q^2/4) from wikipedia
  Real r2 = square(r);
  // this is (-p^3/27) from wikipedia
  Real q3 = cube(q);
//this is the (opposite of the) Cardano assumption ((q^2/4) + (p^3/27)) > 0
//note that r2 may be arbitrarily close to q3, so both
//sides of this if statement have to deal with duplicate roots
  if(r2 < q3) {
// https://en.wikipedia.org/wiki/Cubic_function#Three_real_roots
// the three roots can be projected from three equidistant points on a circle
    q = -2. * sqrt(q);
    /* q is the circle radius, so if it is small enough we have a triple root */
    if (fabs(q) < 1e-5) {
      roots[0] = -a;
      mults[0] = 3;
      return 1;
    }
    Real t = r / sqrt(q3);
    /* ensure we don't exceed the domain of acos()
       (should only exceed due to roundoff) */
    t = max2(-1., t);
    t = min2( 1., t);
    t = acos(t); // angle of rotation
    a /= 3.;
    roots[0] = q * cos(t / 3.) - a;
    roots[1] = q * cos((t + 2. * PI) / 3.) - a;
    roots[2] = q * cos((t - 2. * PI) / 3.) - a;
    mults[0] = 1;
    mults[1] = 1;
    mults[2] = 1;
    /* detect double root cases here: this is just rotation
       of an equilateral triangle to angle (t); there are cases when
       two vertices line up along x (edge is vertical) */
    if (are_close(roots[0], roots[1], 5e-6, 5e-6)) {
      roots[0] = roots[2];
      roots[1] = average(roots[0], roots[1]);
      mults[1] = 2;
      return 2;
    }
    if (are_close(roots[1], roots[2], 5e-6, 5e-6)) {
      roots[1] = average(roots[1], roots[2]);
      mults[1] = 2;
      return 2;
    }
    if (are_close(roots[0], roots[2], 5e-6, 5e-6)) {
      roots[0] = average(roots[0], roots[2]);
      swap2(roots[0], roots[1]);
      mults[1] = 2;
      return 2;
    }
    return 3;
  } else {
// https://en.wikipedia.org/wiki/Cubic_function#Cardano.27s_method
    Real u3 = -r - sqrt(r2 - q3);
  //std::pow will not accept a negative base (it can't tell
  //that (1./3.) is exactly the reciprocal of an odd number),
  //so we could strip out the sign on input and put it back
  //on output:
  //Real u = sign(u3)*pow(fabs(u3), 1./3.);
  //even better, C++11 provides std::cbrt which solves this
    Real u = cbrt(u3);
    Real v = (u == 0.0) ? 0.0 : (q / u);
    Real t1 = u + v;
    // recall x = t - (b/(3a)), in our case x = t - (a/3)
    roots[0] = t1 - (a / 3.);
    Real t_real = -0.5 * (u + v);
    Real t_imag = 0.5 * sqrt(3.) * (u - v);
    roots[1] = (t_real) - (a / 3.);
    if (fabs(t_imag) < 1e-6) {
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

// solve quadratic equation x^2 + a * x + b = 0
INLINE UInt solve_quadratic(
    Real a, Real b,
    Few<Real, 2>& roots, Few<UInt, 2>& mults) {
  Real disc = square(a) - 4. * b;
  if (fabs(disc) < 1e-6) {
    mults[0] = 2;
    roots[0] = -a / 2.;
    return 1;
  }
  if (disc > 0.0) {
    mults[0] = 1;
    mults[1] = 1;
    roots[0] = (-a + sqrt(disc)) / 2.;
    roots[1] = (-a - sqrt(disc)) / 2.;
    return 2;
  }
  return 0;
}
