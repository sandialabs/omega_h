INLINE Real cubic_discriminant(
    Real a, Real b, Real c, Real d) {
  return 18 * a * b * c * d
       -  4 * b * b * b * d
       +      b * b * c * c
       -  4 * a * c * c * c
       - 27 * a * a * d * d;
}

INLINE UInt solve_cubic(
    Real a, Real b, Real c, Real d,
    Few<Real, 3>& roots) {
  Real disc = cubic_discriminant(a,b,c,d);
  UInt nroots = (disc > 0.0) ? 3 : 1;
  Complex u[3];
  u[0] = 1.;
  u[1] = Complex(-1.,  sqrt(3.)) / 2.;
  u[2] = Complex(-1., -sqrt(3.)) / 2.;
  Real d0 = b * b - 3 * a * c;
  Real d1 = 2 * b * b * b
         -  9 * a * b * c
         + 27 * a * a * d;
  bool have_d0 = (fabs(d0) > EPSILON);
  Complex tmp = (have_d0) ? (sqrt(square(d1) - 4 * cube(d0))) : (d1);
  Complex C = pow((d1 + tmp) / 2., 1. / 3.);
  Complex x[3];
  if (abs(C) > EPSILON) {
    for (UInt k = 0; k < 3; ++k) {
      x[k] = -(1. / (3. * a)) *
        (b + u[k] * C + (d0 / (u[k] * C)));
    }
  } else {
    if (have_d0) {
      x[0] = x[1] = (9. * a * d - b * c) / (2. * d0);
      x[2] = (4. * a * b * c - 9. * a * a * d - b * b * b) / (a * d0);
    } else {
      for (UInt k = 0; k < 3; ++k) {
        x[k] = -(b / (3. * a));
      }
    }
  }
  if (nroots == 3) {
    for (UInt k = 0; k < 3; ++k) {
      roots[k] = x[k].real();
    }
  } else {
    UInt best = 0;
    /* return the least imaginary root */
    for (UInt k = 1; k < 3; ++k) {
      if (fabs(x[k].imag()) < fabs(x[best].imag())) {
        best = k;
      }
    }
    roots[0] = x[best].real();
  }
  return nroots;
}
