template <UInt dim>
Real metric_product(
    Matrix<dim,dim> m,
    Vector<dim> v) {
  return v * (m * v);
}

template <UInt dim>
Matrix<dim,dim> compose_metric(
    Matrix<dim,dim> r,
    Vector<3> eigenw) {
  Matrix<dim,dim> l = diagonal(eigenw);
  return r * l * transpose(r);
}

template <UInt dim>
void decompose_metric(
    Matrix<dim,dim> m,
    Matrix<dim,dim>& r,
    Vector<3>& eigenw) {
  Matrix<dim,dim> l;
  decompose_eigen_qr(m, r, l);
  eigenw = diagonal(l);
}

/* INRIA knows what they're doing with respect to the metric.
   See these papers:

   Frey, Pascal-Jean, and Frédéric Alauzet.
   "Anisotropic mesh adaptation for CFD computations."
   Computer methods in applied mechanics and engineering
   194.48 (2005): 5068-5082.

   F. Alauzet, P.J. Frey, Estimateur d'erreur geometrique
   et metriques anisotropes pour l'adaptation de maillage.
   Partie I: aspects theoriques,
   RR-4759, INRIA Rocquencourt, 2003.

https://www.rocq.inria.fr/gamma/Frederic.Alauzet/
*/

template <UInt dim>
Matrix<dim,dim> common_basis(
    Matrix<dim,dim> a,
    Matrix<dim,dim> b) {
  auto c = invert(a) * b;
  Matrix<dim,dim> r;
  Few<Real, 3> eigenw;
  decompose_metric(c, r, eigenw);
  (void) eigenw;
  return r;
}

template <UInt dim>
Matrix<dim,dim> intersect_metric(
    Matrix<dim,dim> a,
    Matrix<dim,dim> b,
    Real t) {
  Matrix<dim,dim> p = common_basis(a, b);
  Few<Real, dim> w;
  for (UInt i = 0; i < dim; ++i) {
    Real u = metric_product(a, p[i]);
    Real v = metric_product(b, p[i]);
    Real w[i] = max2(u, v);
  }
}
