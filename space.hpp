INLINE Vector<2> vector_2(Real x, Real y) {
  Vector<2> v;
  v[0] = x; v[1] = y;
  return v;
}

INLINE Vector<3> vector_3(Real x, Real y, Real z) {
  Vector<3> v;
  v[0] = x; v[1] = y; v[2] = z;
  return v;
}

INLINE Matrix<2,2> matrix_2x2(
    Real a, Real b,
    Real c, Real d) {
  Matrix<2,2> o;
  o[0] = vector_2(a, c);
  o[1] = vector_2(b, d);
  return o;
}

INLINE Matrix<3,3> matrix_3x3(
    Real a, Real b, Real c,
    Real d, Real e, Real f,
    Real g, Real h, Real i) {
  Matrix<3,3> o;
  o[0] = vector_3(a, d, g);
  o[1] = vector_3(b, e, h);
  o[2] = vector_3(c, f, i);
  return o;
}

INLINE Matrix<3,3> cross(Vector<3> a) {
  return matrix_3x3(
        0  , -a[2],  a[1],
       a[2],   0  , -a[0],
      -a[1],  a[0],   0  );
}

INLINE Vector<3> cross(Vector<3> a, Vector<3> b) {
  return vector_3(
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]);
}

INLINE Real cross(Vector<2> a, Vector<2> b) {
  return (a[0] * b[1] - a[1] * b[0]);
}

/* Rodrigues' Rotation Formula */
INLINE Matrix<3,3> rotate(Real angle, Vector<3> axis) {
  return cos(angle) * identity_matrix<3,3>() +
         sin(angle) * cross(axis) +
         (1 - cos(angle)) * tensor_product(axis, axis);
}

INLINE Matrix<2,2> rotate(Real angle) {
  return matrix_2x2(
      cos(angle), -sin(angle),
      sin(angle),  cos(angle));
}

INLINE Real det(Matrix<2,2> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return a * d - b * c;
}

INLINE Real det(Matrix<3,3> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[2][0];
  Real d = m[0][1];
  Real e = m[1][1];
  Real f = m[2][1];
  Real g = m[0][2];
  Real h = m[1][2];
  Real i = m[2][2];
  return (a * e * i)
       + (b * f * g)
       + (c * d * h)
       - (c * e * g)
       - (b * d * i)
       - (a * f * h);
}

INLINE Matrix<2,2> invert(Matrix<2,2> m) {
  Real a = m[0][0];
  Real b = m[1][0];
  Real c = m[0][1];
  Real d = m[1][1];
  return matrix_2x2(
      d, -b,
     -c,  a) / det(m);
}

INLINE Matrix<3,3> invert(Matrix<3,3> a) {
  Matrix<3,3> b;
  b[0] = cross(a[1], a[2]);
  b[1] = cross(a[2], a[0]);
  b[2] = cross(a[0], a[1]);
  return transpose(b) / det(a);
}
