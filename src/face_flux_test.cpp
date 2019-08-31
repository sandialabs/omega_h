#include <Omega_h_build.hpp>
#include <Omega_h_coarsen.hpp>
#include <Omega_h_refine.hpp>
#include <Omega_h_compare.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_metric.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_reduce.hpp>
#include <Omega_h_assoc.hpp>
#include <Omega_h_expr.hpp>
#include <Omega_h_align.hpp>

#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>

using namespace Omega_h;

void FunctionInitialCondition(
    const MeshDimSets &facesets,
    const MeshDimSets &elemsets,
    const Reals       &coords,
    const LOs         &face_node_id,
    const LOs         &elem_face_id,
    Write<Real>        field);


void FunctionInitialCondition(
    const MeshDimSets &facesets,
    const MeshDimSets &elemsets,
    const Reals       &coords,
    const LOs         &face_node_id,
    const LOs         &elem_face_id,
    Write<Real>        field) {

  constexpr int N = 3;
  constexpr int F = 4;
  constexpr int D = 3;
  const std::vector<std::string> 
    face_sets={"x-","y-","z-","x+","y+","z+"};
  const std::vector<std::string>
    element_blocks={"body"};
  const std::string expr_string="vector(2*z, 5*x, y)";
  const std::string name="Face flux test";

  for (const std::string &faceset : face_sets) {
    auto fsIter = facesets.find(faceset);
    if(fsIter == facesets.end())
       fail("Faceset block doesn't exist!");
    const auto     faceLids = (fsIter->second);
    const auto   nset_faces = faceLids.size();
    // Edge Centers
    Write<double> Cx(N*nset_faces);
    Write<double> Cy(N*nset_faces);
    Write<double> Cz(N*nset_faces);
    // Edge Vectors
    Write<double> Vx(N*nset_faces);
    Write<double> Vy(N*nset_faces);
    Write<double> Vz(N*nset_faces);

    const Reals cord = coords;
    const LOs    ids = face_node_id;
    auto prepare = OMEGA_H_LAMBDA(int set_face) {
      const auto face = faceLids[set_face];
      int nodes[N];
      for (int i = 0; i < N; ++i)
        nodes[i] = ids[N*face+i];
      double X[N][3];
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
           X[i][j] = cord[D*nodes[i]+j];

      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i) {
        const int j = (i+1)%N;
        Vx[offset+i] =  X[j][0] - X[i][0];
        Vy[offset+i] =  X[j][1] - X[i][1];
        Vz[offset+i] =  X[j][2] - X[i][2];
        Cx[offset+i] = (X[j][0] + X[i][0])/2;
        Cy[offset+i] = (X[j][1] + X[i][1])/2;
        Cz[offset+i] = (X[j][2] + X[i][2])/2;
      }
    };
    parallel_for(nset_faces, prepare);
    const std::vector<std::string> S={"x","y","z"};
    ExprReader reader(N*nset_faces, D);
    reader.register_variable(S[0], Omega_h::any(Reals(Cx)));
    reader.register_variable(S[1], Omega_h::any(Reals(Cy)));
    reader.register_variable(S[2], Omega_h::any(Reals(Cz)));
    auto result = reader.read_string(expr_string, name);
    reader.repeat(result);
    Reals field_osh = any_cast<Reals>(result);

    Write<Real> fld=field;
    auto save = OMEGA_H_LAMBDA(int set_face) {
      double Fl[N][3];
      double v [N][3];
      const auto face = faceLids[set_face];
      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i)
        for (int s = 0; s < D; ++s)
          Fl[i][s] = field_osh[(offset+i) * D + s];
      for (int i = 0; i < N; ++i) {
        v[i][0] = Vx[offset+i];
        v[i][1] = Vy[offset+i];
        v[i][2] = Vz[offset+i];
      }
      double flux=0;
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
          flux += Fl[i][j]*v[i][j];
      fld[face] = flux;
    };
    parallel_for(nset_faces, save);
  }

  for (const std::string &blockname : element_blocks) {
    auto esIter = elemsets.find(blockname);
    if(esIter == elemsets.end())
       fail("Element block doesn't exist!");
    const auto elementLids = (esIter->second);
    const int   nset_elems = elementLids.size();
    const int  nset_faces = F*nset_elems;
    // Edge Centers
    Write<double> Cx(N*nset_faces);
    Write<double> Cy(N*nset_faces);
    Write<double> Cz(N*nset_faces);
    // Edge Vectors
    Write<double> Vx(N*nset_faces);
    Write<double> Vy(N*nset_faces);
    Write<double> Vz(N*nset_faces);
    const Reals cord = coords;
    const LOs   ids  = face_node_id;
    const LOs   fids = elem_face_id;
    auto prepare = OMEGA_H_LAMBDA(int set_face) {
      const int  el = elementLids[set_face/F];
      const int   f = set_face%F;
      const auto face = fids[F*el+f];
      int nodes[N];
      for (int i = 0; i < N; ++i)
        nodes[i] = ids[N*face+i];
      double X[N][3];
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
           X[i][j] = cord[D*nodes[i]+j];

      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i) {
        const int j = (i+1)%N;
        Vx[offset+i] =  X[j][0] - X[i][0];
        Vy[offset+i] =  X[j][1] - X[i][1];
        Vz[offset+i] =  X[j][2] - X[i][2];
        Cx[offset+i] = (X[j][0] + X[i][0])/2;
        Cy[offset+i] = (X[j][1] + X[i][1])/2;
        Cz[offset+i] = (X[j][2] + X[i][2])/2;
      }
    };
    parallel_for(nset_faces, prepare);
    const std::vector<std::string> S={"x","y","z"};
    ExprReader reader(N*nset_faces, D);
    reader.register_variable(S[0], Omega_h::any(Reals(Cx)));
    reader.register_variable(S[1], Omega_h::any(Reals(Cy)));
    reader.register_variable(S[2], Omega_h::any(Reals(Cz)));
    auto result = reader.read_string(expr_string, name);
    reader.repeat(result);
    Reals field_osh = any_cast<Reals>(result);

    Write<Real> fld=field;
    auto save = OMEGA_H_LAMBDA(int set_face) {
      double Fl[N][3];
      double v [N][3];
      const int  el = elementLids[set_face/F];
      const int   f = set_face%F;
      const auto face = fids[F*el+f];
      const auto offset = set_face*N;
      for (int i = 0; i < N; ++i)
        for (int s = 0; s < D; ++s)
          Fl[i][s] = field_osh[(offset+i) * D + s];
      for (int i = 0; i < N; ++i) {
        v[i][0] = Vx[offset+i];
        v[i][1] = Vy[offset+i];
        v[i][2] = Vz[offset+i];
      }
      double flux=0;
      for (int i = 0; i < N; ++i)
        for (int j = 0; j < D; ++j)
          flux += Fl[i][j]*v[i][j];

      fld[face] = flux;
    };
    parallel_for(nset_faces, save);
  }
}

OMEGA_H_INLINE
Few<Vector<3>,4>
comp_reference_nodal_gradients(const double *const /*xi*/)
{
  typedef Vector<3> Vector;
  typedef Few<Vector,4> Bucket;
  Bucket returnMe;
  {
    Vector &GRAD = returnMe[0];
    GRAD(0) = -1.;
    GRAD(1) = -1.;
    GRAD(2) = -1.;
  }
  {
    Vector &GRAD = returnMe[1];
    GRAD(0) = +1.;
    GRAD(1) = +0.;
    GRAD(2) = +0.;
  }
  {
    Vector &GRAD = returnMe[2];
    GRAD(0) = +0.;
    GRAD(1) = +1.;
    GRAD(2) = +0.;
  }
  {
    Vector &GRAD = returnMe[3];
    GRAD(0) = +0.;
    GRAD(1) = +0.;
    GRAD(2) = +1.;
  }
  return returnMe;
}

OMEGA_H_INLINE
Few<Vector<3>,4>
comp_face_basis( const double *const x,
                 const double *const y,
                 const double *const z,
                 const double *const xi )
{
  typedef Vector<3> Vector;
  typedef Matrix<3,3> Tensor;
  typedef Few<Vector,4> Bucket;

  const Bucket gradN = comp_reference_nodal_gradients(xi);
  Tensor JinvF = zero_matrix<3,3>();
  for (int A=0; A<4; ++A) {
    Vector xA;
    xA(0) = x[A];
    xA(1) = y[A];
    xA(2) = z[A];
    Vector gradA = gradN[A];
    const Tensor addMe = outer_product(xA,gradA);
    JinvF += addMe;
  }
  const double J = determinant(JinvF);
  JinvF /= J;

  Bucket returnMe;
  {
    //Omega_h face 0 = intrepid face +3
    Vector &fillMe = returnMe[0];
    fillMe(0) = 2.0*xi[0];
    fillMe(1) = 2.0*xi[1];
    fillMe(2) = 2.0*(xi[2] - 1.0);
    fillMe = JinvF*fillMe;
  }
  {
    //Omega_h face 1 = intrepid face +0
    Vector &fillMe = returnMe[1];
    fillMe(0) = 2.0*xi[0];
    fillMe(1) = 2.0*(xi[1] - 1.0);
    fillMe(2) = 2.0*xi[2];
    fillMe = JinvF*fillMe;
  }
  {
    //Omega_h face 2 = intrepid face +1
    Vector &fillMe = returnMe[2];
    fillMe(0) = 2.0*xi[0];
    fillMe(1) = 2.0*xi[1];
    fillMe(2) = 2.0*xi[2];
    fillMe = JinvF*fillMe;
  }
  {
    //Omega_h face 3 = intrepid face +2
    Vector &fillMe = returnMe[3];
    fillMe(0) = 2.0*(xi[0] - 1.0);
    fillMe(1) = 2.0*xi[1];
    fillMe(2) = 2.0*xi[2];
    fillMe = JinvF*fillMe;
  }

  return returnMe;
}

namespace {
typedef std::array<double,3> V;
double dot(const V &a, const V &b) {
  double r = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  return r;
}
V cross(const V &a, const V &b) {
  V v={{a[1]*b[2] - a[2]*b[1],
       -a[0]*b[2] + a[2]*b[0],
        a[0]*b[1] - a[1]*b[0]}};
  return v;
}
double norm(const V &v) {
  const double r = dot(v,v);
  return std::sqrt(r); 
}
}
int main(int argc, char** argv) {
  Library lib(&argc, &argv);
  CommPtr world = lib.world();
  Mesh mesh = build_box(world, OMEGA_H_SIMPLEX, 1., 2., 3., 2, 2, 2);
  AdaptOpts opts(&mesh);

  mesh.add_tag<Real>(VERT, "metric", 1);
  const std::string name = "magnetic_face_flux";

  {
    Write<Real> OK(1.,0);
    Assoc assoc = get_box_assoc(3);
    const MeshSets mesh_sets = invert(&mesh, assoc);
    const MeshDimSets facesets = mesh_sets[SIDE_SET];
    const MeshDimSets elemsets = mesh_sets[ELEM_SET];
    LOs face_node_id = mesh.ask_down(FACE, VERT).ab2b;
    LOs elem_node_id = mesh.ask_down(REGION, VERT).ab2b;
    auto const ElemFace = mesh.ask_down(3,2);
    LOs elem_face_id = ElemFace.ab2b;

    Reals coords  = mesh.coords();
    Write<Real> flux_w(mesh.nfaces());
    FunctionInitialCondition(
      facesets,
      elemsets,
      coords,
      face_node_id,
      elem_face_id,
      flux_w);

    Reals flux(flux_w);
    mesh.add_tag<Real>(FACE, name, 1, flux);

    const std::vector<std::string> element_blocks={"body"};
    for (const std::string &blockname : element_blocks) {
      constexpr int N = 4;
      constexpr int F = 4;
      constexpr int D = 3;
      auto esIter = elemsets.find(blockname);
      if(esIter == elemsets.end())
         fail("Element block doesn't exist!");
      const int   nset_elems = esIter->second.size();
      const Reals cord = coords;
      const LOs   nids  = elem_node_id;
      const LOs   fids  = elem_face_id;
      Read<Real> flux_r=mesh.get_array<Real>(FACE, name); 
      auto check = OMEGA_H_LAMBDA(int elem) {
        int nodes[N];
        int faces[F];
        for (int i = 0; i < N; ++i) nodes[i] = nids[N*elem+i];
        for (int i = 0; i < F; ++i) faces[i] = fids[F*elem+i];
        double X[D][N];
        for (int i = 0; i < N; ++i)
          for (int j = 0; j < D; ++j)
             X[j][i] = cord[D*nodes[i]+j];
        constexpr double xi[] = {1./3.,1./3.,1./3.};
        const Few<Vector<3>,4> faceBasis =
          comp_face_basis(X[0],X[1],X[2],xi);
        Vector<3> B = zero_vector<3>();
        for (int f=0; f<F; ++f) {
          const I8 code = ElemFace.codes[elem*F+f];
          const int sign = code_is_flipped(code) ? -1 : +1;
          B += sign * flux_r[faces[f]] * faceBasis[f];
        }
        const Vector<3> y = {{1, 2, 5}};
        const double tol = 1.0e-12;
        for (int i = 0; i < D; ++i)
          if (tol < std::abs(B[i] - y[i])) 
            OK[0] += 1;
      };
      parallel_for(nset_elems, check);
    }
    const bool ok = 0.==Reals(OK)[0];
    if (!ok) return 2;
  }



  mesh.set_tag(
      VERT, "metric", Reals(mesh.nverts(), metric_eigenvalue_from_length(1.3)));
  while (coarsen_by_size(&mesh, opts))
    ;



  {
    Write<Real> OK(1.,0);
    Assoc assoc = get_box_assoc(3);
    const MeshSets mesh_sets = invert(&mesh, assoc);
    const MeshDimSets elemsets = mesh_sets[ELEM_SET];
    LOs elem_node_id = mesh.ask_down(REGION, VERT).ab2b;
    auto const ElemFace = mesh.ask_down(3,2);
    LOs elem_face_id = ElemFace.ab2b;
    Reals coords  = mesh.coords();
  
    const std::vector<std::string> element_blocks={"body"};
    for (const std::string &blockname : element_blocks) {
      constexpr int N = 4;
      constexpr int F = 4;
      constexpr int D = 3;
      auto esIter = elemsets.find(blockname);
      if(esIter == elemsets.end())
         fail("Element block doesn't exist!");
      const int   nset_elems = esIter->second.size();
      const Reals cord = coords;
      const LOs   nids  = elem_node_id;
      const LOs   fids  = elem_face_id;
      Read<Real> flux_r=mesh.get_array<Real>(FACE, name); 
      auto check = OMEGA_H_LAMBDA(int elem) {
        int nodes[N];
        int faces[F];
        for (int i = 0; i < N; ++i) nodes[i] = nids[N*elem+i];
        for (int i = 0; i < F; ++i) faces[i] = fids[F*elem+i];
        double X[D][N];
        for (int i = 0; i < N; ++i)
          for (int j = 0; j < D; ++j)
             X[j][i] = cord[D*nodes[i]+j];
        constexpr double xi[] = {1./3.,1./3.,1./3.};
        const Few<Vector<3>,4> faceBasis =
          comp_face_basis(X[0],X[1],X[2],xi);
        Vector<3> B = zero_vector<3>();
        for (int f=0; f<F; ++f) {
          const I8 code = ElemFace.codes[elem*F+f];
          const int sign = code_is_flipped(code) ? -1 : +1;
          B += sign * flux_r[faces[f]] * faceBasis[f];
        }
        const Vector<3> y = {{1, 2, 5}};
        const double tol = 1.0e-12;
        for (int i = 0; i < D; ++i)
          if (tol < std::abs(B[i] - y[i])) 
            OK[0] += 1;

        double tot_f = 0;
        for (int f = 0; f < F; ++f) {
          typedef std::array<double,D> V;
          std::array<V,3> x;
          for (int i = 0,n=0; i < N; ++i)
            if (i != f) {
              for (int j = 0; j < D; ++j)
                x[n][j] = X[j][i];
              ++n;
            }
          const V v1 = {{x[1][0]-x[0][0],x[1][1]-x[0][1],x[1][2]-x[0][2]}};
          const V v2 = {{x[2][0]-x[0][0],x[2][1]-x[0][1],x[2][2]-x[0][2]}};
          const V v3 = {{X[0][f]-x[0][0],X[1][f]-x[0][1],X[2][f]-x[0][2]}};
          V n = cross(v1,v2);
          const double sign = -std::copysign(1.0,dot(n,v3));
          const double a = norm(n);
          n[0] /= sign*a;
          n[1] /= sign*a;
          n[2] /= sign*a;
          const V g = {{1,2,5}};
          const double Ft = dot(n,g) * a/2;
          const int h = f==0 ? 2 : f==1 ? 3 : f==2 ? 1 : 0;
          const I8 code = ElemFace.codes[elem*F+h];
          const int f_sign = code_is_flipped(code) ? -1 : +1;
          const double ft = f_sign*flux_r[faces[h]];
          tot_f += ft;
          if (tol < std::abs(Ft-ft))  OK[0] += 1;
        }
        if (.000001 < std::abs(tot_f) ) OK[0] += 1;
      };
      parallel_for(nset_elems, check);
    }
    const bool ok = 0.==Reals(OK)[0];
    if (!ok) return 2;
  }
  return 0;
}
